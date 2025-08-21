# -*- coding: utf-8 -*-
"""
Created on Thu May 29 10:04:24 2025

@author: hgtra
"""

from astropy.io import fits
from astropy.table import Table, vstack
import os
import numpy as np
from astropy.coordinates import SkyCoord, match_coordinates_sky
import time
from tqdm import tqdm
import gc

# can be set to 1 if MemoryError, to speed up subsequent iterations.
exclude_distant = 0
# will discard new matches to D>2500 galaxies 
# => decrease distance reliability, less ellipses

# Catalog to start with. If >0, resume last run
start = 0
if start>9:
    exclude_distant=1
no_id = 0

home = os.path.expanduser("~")
os.chdir(home+"/Downloads/Catalogs/regalade")

catalogs  =  ['sga_input.fits', #0, 1
              'glade1_input.fits', #1, 2
              'hecate_input.fits', #2, 4
              'desi_pv_input.fits', #3, 8
              'desi_dr1_input.fits', #4, 16
              'cosmicflows_input.fits', #5, 32
              'nedlvs_d_input.fits', #6, 64 
              'nedlvs_z_input.fits', #7, 128
              'nedlvs_r_input.fits', #8, 256
              'gladep_input.fits', #9, 512 (1st catalog with >1e7 rows)
              'desi_zou_input.fits', #10, 1024
              'panstarrs_input.fits', #11, 2048
              'sdssgal_input.fits', #12, 4096
              'gsc_input.fits', #13, 8192
              'ls10_wh24_input.fits', #14, 16384
              'delve_input.fits', #15, 32768
              ]

# start with small catalogs, ideally with ellipses, should be all kept
# before the compilation exceeds 1e7 rows, keep all rows: some just have missing distance
# when the compilation exceeds 1e7 rows, only D>2000 orphans are added
# priority to catalogs with ellipses and reliable D / large fraction not null
# use gladep / ls10 before delve => delve morph added even if null D
# finish with deep surveys
# order dist:  5, 6, 7, 4, 8, 0, 2, 3, 9, 1, 11, 13, 10, 14, 15, 12,
# (used in next code compute_dist.py)

tmatch_cols = [f"t{j}_match_idx" for j in range(len(catalogs))]

def match_catalogs(i, table_a, table_b, c_a, c_b):
    print(len(table_a),len(table_b))

    if len(table_a) < 1e7 or len(table_a) < len(table_b):
        # Match a to b directly
        idx, sep2d, _ = match_coordinates_sky(c_a, c_b)
        sep2d = sep2d.arcsec
    else:
        print("reverse match")
        # Reverse match: b to a
        ind, d2d, _ = match_coordinates_sky(c_b, c_a)
        d2d = d2d.arcsec
        
        # Step 2: Sort by distance (closest first)
        isort = np.argsort(d2d)
        ind_sorted = ind[isort]
        d2d_sorted = d2d[isort]

        # Step 3: Find unique a matches (closest b for each a)
        _, iu = np.unique(ind_sorted, return_index=True)
        keep_mask = np.zeros(len(ind_sorted), dtype=bool)
        keep_mask[iu] = True

        # Step 4: Reverse the sorting to original b order
        rsort = np.argsort(isort)
        keep_mask = keep_mask[rsort]

        # Step 5: Build final index array (length = len(c_a))
        idx = np.full(len(c_a), -1, dtype=int)
        sep2d = np.full(len(c_a), np.inf)

        # For all valid b[i] with unique closest a, assign
        a_indices = ind[keep_mask]
        b_indices = np.nonzero(keep_mask)[0]
        
        idx[a_indices] = b_indices
        sep2d[a_indices] = d2d[keep_mask]
        
        
        
    matched = sep2d < 5
    
    # Compute position angle difference
    posang_a = c_b[idx].position_angle(c_a).deg - np.asarray(table_a['PA'])
    posang_b = c_a.position_angle(c_b[idx]).deg - np.asarray(table_b[idx]['PA'])

    # Apply elliptical containment condition
    R1_a = np.asarray(table_a['R1'])
    R2_a = np.asarray(table_a['R2'])
    R1_b = np.asarray(table_b[idx]['R1'])
    R2_b = np.asarray(table_b[idx]['R2'])
    reassign_ellipse = np.logical_and(matched,
                                      np.logical_or(
                                          np.logical_and(i==15, np.logical_or(R1_a<R1_b,R1_b>6)), # delve
                                          np.logical_and(R1_a==3,R1_b>3) # new morph
                                          ))
    table_a['ref_ellipse'][reassign_ellipse] = i
    
    # Redefine galaxy ellipse (a better morpho was found)
    table_a['R1'][reassign_ellipse] = table_b[idx[reassign_ellipse]]['R1']
    table_a['R2'][reassign_ellipse] = table_b[idx[reassign_ellipse]]['R2']
    table_a['PA'][reassign_ellipse] = table_b[idx[reassign_ellipse]]['PA']
    table_a['ra'][reassign_ellipse] = table_b[idx[reassign_ellipse]]['ra']
    table_a['dec'][reassign_ellipse] = table_b[idx[reassign_ellipse]]['dec']
    
    
    factor = 1
    rad_min = np.maximum(factor * np.maximum(
                                            np.minimum((R1_a * R2_a / np.sqrt(R1_a**2 * np.sin(np.radians(posang_a))**2 + R2_a**2 * np.cos(np.radians(posang_a))**2)),
                                                       (R1_b * R2_b / np.sqrt(R1_b**2 * np.sin(np.radians(posang_b))**2 + R2_b**2 * np.cos(np.radians(posang_b))**2))),
                                            (1/4.)*(R1_a * R2_a / np.sqrt(R1_a**2 * np.sin(np.radians(posang_a))**2 + R2_a**2 * np.cos(np.radians(posang_a))**2))),
                         5.)
    matched = sep2d < rad_min
    if np.any(idx[matched]<0):
        raise SystemExit("fatal error: reverse match failed")

    return table_a, matched, table_b['t%d_match_idx'%i][idx[matched]]



for i in tqdm(range(start,len(catalogs))):
    print(catalogs[i])
    if 'gladep' in catalogs[i-1]:
        exclude_distant = 1
    if i == start:
        # Initialization
        if start == 0:
            master = Table.read(catalogs[i])
            master['ref_ellipse'] = 0
            master['t0_match_idx'] = np.arange(len(master))+1
            
        else:
            # force resuming last run
            master = Table.read("master_catalog%s.fits"%("_noID"*no_id))
            
        coo0 = SkyCoord(master['ra'],master['dec'],frame='icrs',unit='deg')
        continue
    
    if not(np.all(master['ra'][:-1] <= master['ra'][1:])):
        ind = np.argsort(master['ra'])
        master = master[ind]
        coo0 = coo0[ind]
        
    # Load secondary catalog
    aux = Table.read(catalogs[i])
    
    if not(np.all(aux['ra'][:-1] <= aux['ra'][1:])):
        raise SystemExit(catalogs[i],"is not sorted")
        aux = aux[np.argsort(aux['ra'])]
    aux['t%d_match_idx'%i] = np.arange(len(aux))+1
    if exclude_distant:
        if np.ma.getmask(aux['D']).any():
            ind = np.logical_or(aux['D'].mask,
                                aux['D']<2500)
        else:
            ind = aux['D']<2500
            
        aux = aux[ind]
        print(len(ind), "->", len(aux))
    
    coo1 = SkyCoord(aux['ra'],aux['dec'],frame='icrs',unit='deg')
    t0 = time.time()
    
    # To avoid MemoryError, split into chunks: e.g. ra<180 and ra>=180
    # All catalogs should be sorted by RA. No check is made, though
    match = np.full(len(master),False)
    nchunks = 4
    parts = []
    idx = []
    
    for k in range(nchunks):
        if k < nchunks-1:
            ind0 = np.logical_and(master['ra']>=360*k/nchunks,
                                  master['ra']<360*(k+1)/nchunks)
            ind1 = np.logical_and(aux['ra']>=360*k/nchunks,
                                  aux['ra']<360*(k+1)/nchunks)
        else:
            # include 360
            ind0 = np.logical_and(master['ra']>=360*k/nchunks,
                                  master['ra']<=360*(k+1)/nchunks)
            ind1 = np.logical_and(aux['ra']>=360*k/nchunks,
                                  aux['ra']<=360*(k+1)/nchunks)
            
        
            
        table_0 = master[ind0]
        table_1 = aux[ind1]
        
        table_0, matched, idx_part = match_catalogs(i, table_0, table_1, 
                                                        coo0[ind0], coo1[ind1])
        parts.append(table_0)
        idx = np.concatenate([idx,idx_part])
        
            
        match[ind0] = matched
        
        t1 = time.time()    
        print("time to match table %d (%d/%d): %.1f"%(i, k+1, nchunks, t1-t0))
        t0 = time.time()
    
    master = vstack(parts)
   
    master['t%d_match_idx'%i] = 0
    master['t%d_match_idx'%i][match] = idx # 1-based
    nonmatch = ~np.isin(aux['t%d_match_idx'%i],idx)
    if len(aux)>1e7:    
        keep = np.logical_and(nonmatch,
                              aux['D']<2000)
    else:
        keep = nonmatch
    
    aux_keep = aux[keep]
    aux_keep['ref_ellipse'] = i    
    for c in master.colnames:
        if not(c in aux_keep.colnames):
            aux_keep[c] = 0
            
    master = vstack([master,aux_keep])
    t1 = time.time()
    print("\nstacking with aux: %.1f seconds"%(t1-t0))
    
    if i==len(catalogs)-1:
        # last iteration: add distance columns
        coo0 = SkyCoord(master['ra'],master['dec'],frame='icrs',unit='deg')
        master['D'].name = 'D_input'
        master['ref_D_input'] = 0
        
        ncat = np.full(len(master),0,'int')
        for c in tmatch_cols:
            ncat[master[c]>0]+=1
        
        # for large catalogs, exclude distant orphans
        zou_exclude = np.all((master['ref_ellipse']==10,
                              ncat==1,
                              np.abs(master['D_input']-700)>500), axis=0)
        
        gsc_exclude = np.all((master['ref_ellipse']==13,
                              ncat==1,
                              np.logical_or(np.abs(coo0.galactic.b.deg)<20,
                                            master['D_input']<200)), axis=0)
        
        wh24_exclude = np.all((master['ref_ellipse']==14,
                              ncat==1,
                              np.abs(master['D_input']-700)>500), axis=0)
        
        exclude = np.any((zou_exclude,
                          gsc_exclude,
                          wh24_exclude),axis=0)
        
        master = master[~exclude]
        print(len(exclude),'->',len(master))
        
        master['D_tmean'] = np.full(len(master), np.nan, 'float')
        master['n_dist'] = 0
        master['id_cat'] = ncat[~exclude]
        
        
    # SAVE TABLE FUNCTION
    cols = []
    for name in master.colnames:
        if no_id and name in tmatch_cols:
            continue
        print(name)
        data = master[name]
        if "idx" in name and np.ma.getmask(data).any():
            data.mask = False
        dtype = data.dtype

        if np.issubdtype(dtype, np.floating):
            format_code = 'E'  # 32-bit float
        elif np.issubdtype(dtype, np.integer):
            format_code = 'J'  # 32-bit integer
        else:
            raise TypeError(f"Unsupported column type for '{name}': {dtype}")

        col = fits.Column(name=name, array=data, format=format_code)
        cols.append(col)

    hdu = fits.BinTableHDU.from_columns(cols)
    del aux
    del aux_keep
    del coo0
    del coo1
    del table_0
    del table_1
    gc.collect()

    hdu.writeto('test.fits', overwrite=True)
    del hdu
    
    t0 = time.time()
    print("saved in %.1f seconds"%(t0-t1))
    os.replace("test.fits","master_catalog%s.fits"%("_noID"*no_id))
    
    coo0 = SkyCoord(master['ra'],master['dec'],frame='icrs',unit='deg')
    