# -*- coding: utf-8 -*-
"""
Created on Thu Apr 24 10:31:48 2025

@author: hgtra
"""

from astropy.io import fits
import os
import numpy as np
from astropy.coordinates import SkyCoord, match_coordinates_sky
import time
t0 = time.time()

# takes about 35 minutes to run

home = os.path.expanduser("~")
os.chdir(home+"/Downloads/Catalogs/regalade")


fname = "master_catalog_thin_without_duplicates.fits"
cat_out = 'master_catalog_thin_enhanced.fits'

with fits.open(fname, memmap=True) as hdul:
    t = hdul[1].data
   
    
smallcat = len(t)<1e4
print("small catalog:",smallcat)

coo = SkyCoord(t['ra'].astype('float32'),t['dec'].astype('float32'),frame='icrs',unit='deg')
    
cols = {'G': np.full(len(t),np.nan,'float'),
        'BP': np.full(len(t),np.nan,'float'),
        'PM': np.full(len(t),np.nan,'float'),
        'angDist': np.full(len(t),np.nan,'float'),
        'rmag': np.full(len(t),np.nan,'float'),
        'gKmag': np.full(len(t),np.nan,'float'),
        'rKmag': np.full(len(t),np.nan,'float'),
        'iKmag': np.full(len(t),np.nan,'float'),
        'zKmag': np.full(len(t),np.nan,'float'),
        'W1mag': np.full(len(t),np.nan,'float'),
        'W2mag': np.full(len(t),np.nan,'float'),
        'dK': np.full(len(t),np.nan,'float'),
        'griz_ref': np.full(len(t),0,'int8'),
        'W_ref': np.full(len(t),0,'int8'),
        'contam_type': np.full(len(t),0,'int8'),
        }
    


### MATCH TO CDS Gaia
gaia = "stacked_cds_gaia.fits"
with fits.open(gaia, memmap=True) as hdul:
    g = hdul[1].data

coo_g = SkyCoord(g['ra'],g['dec'],frame='icrs',unit='deg')
if smallcat:
    ind, sep2d, _ = match_coordinates_sky(coo, coo_g, storekdtree=False)
    sep2d = sep2d.arcsec
    # Step 2: Sort by distance (closest first)
    isort = np.argsort(sep2d)
    ind_sorted = ind[isort]
    # Step 3: Find unique a matches (closest b for each a)
    _, iu = np.unique(ind_sorted, return_index=True)
    keep_mask = np.zeros(len(ind_sorted), dtype=bool)
    keep_mask[iu] = True
    # Step 4: Reverse the sorting to original b order
    rsort = np.argsort(isort)
    keep_mask = keep_mask[rsort]
    # Step 5: Build final index array (length = len(c_a))
    idx = np.full(len(coo_g), -1, dtype=int)
    d2d = np.full(len(coo_g), np.inf)
    # For all valid b[i] with unique closest a, assign
    a_indices = ind[keep_mask]
    b_indices = np.nonzero(keep_mask)[0]
    idx[a_indices] = b_indices
    d2d[a_indices] = sep2d[keep_mask]
    match = d2d < 3
else:    
    idx,d2d,_ = coo_g.match_to_catalog_sky(coo)
    match = d2d.arcsec < 3
g = g[match]
for cname in ['Gmag_x','BPmag','PM','angDist']:
    c = cname.replace('mag_x','').replace('mag','')
    cols[c][idx[match]] = g[cname]
    
t1 = time.time()
print("matched to Gaia after %.1f seconds"%(t1-t0))
t0 = time.time()

### MATCH TO CDS PanSTARRS, AllWISE
griz = "stacked_cds_grizW1W2.fits"
with fits.open(griz, memmap=True) as hdul:
    g = hdul[1].data

coo_g = SkyCoord(g['ra'],g['dec'],frame='icrs',unit='deg')
if smallcat:
    ind, sep2d, _ = match_coordinates_sky(coo, coo_g, storekdtree=False)
    sep2d = sep2d.arcsec
    # Step 2: Sort by distance (closest first)
    isort = np.argsort(sep2d)
    ind_sorted = ind[isort]
    # Step 3: Find unique a matches (closest b for each a)
    _, iu = np.unique(ind_sorted, return_index=True)
    keep_mask = np.zeros(len(ind_sorted), dtype=bool)
    keep_mask[iu] = True
    # Step 4: Reverse the sorting to original b order
    rsort = np.argsort(isort)
    keep_mask = keep_mask[rsort]
    # Step 5: Build final index array (length = len(c_a))
    idx = np.full(len(coo_g), -1, dtype=int)
    d2d = np.full(len(coo_g), np.inf)
    # For all valid b[i] with unique closest a, assign
    a_indices = ind[keep_mask]
    b_indices = np.nonzero(keep_mask)[0]
    idx[a_indices] = b_indices
    d2d[a_indices] = sep2d[keep_mask]
    match = d2d < 1
else:
    idx,d2d,_ = coo_g.match_to_catalog_sky(coo)
    match = d2d.arcsec < 1
g = g[match]
for c in cols.keys():
    if c in g.columns.names:
        if c=='W1mag':
            # conversion to AB magnitudes
            cols[c][idx[match]] = g[c]+2.699
        elif c=='W2mag':
            cols[c][idx[match]] = g[c]+3.339
        else:
            cols[c][idx[match]] = g[c]

idx_griz = idx[match][~np.isnan(cols['rmag'][idx[match]])]
idx_w = idx[match][~np.isnan(cols['W1mag'][idx[match]])]
cols['griz_ref'][idx_griz] = 1
cols['W_ref'][idx_w] = 1

t1 = time.time()
print("matched to grizW1W2 after %.1f seconds"%(t1-t0))
t0 = time.time()

### MATCH TO CDS AllWISE, 2MASS XSC
ks = "stacked_cds_k.fits"
with fits.open(ks, memmap=True) as hdul:
    g = hdul[1].data

coo_g = SkyCoord(g['ra'],g['dec'],frame='icrs',unit='deg')
if smallcat:
    ind, sep2d, _ = match_coordinates_sky(coo, coo_g, storekdtree=False)
    sep2d = sep2d.arcsec
    # Step 2: Sort by distance (closest first)
    isort = np.argsort(sep2d)
    ind_sorted = ind[isort]
    # Step 3: Find unique a matches (closest b for each a)
    _, iu = np.unique(ind_sorted, return_index=True)
    keep_mask = np.zeros(len(ind_sorted), dtype=bool)
    keep_mask[iu] = True
    # Step 4: Reverse the sorting to original b order
    rsort = np.argsort(isort)
    keep_mask = keep_mask[rsort]
    # Step 5: Build final index array (length = len(c_a))
    idx = np.full(len(coo_g), -1, dtype=int)
    d2d = np.full(len(coo_g), np.inf)
    # For all valid b[i] with unique closest a, assign
    a_indices = ind[keep_mask]
    b_indices = np.nonzero(keep_mask)[0]
    idx[a_indices] = b_indices
    d2d[a_indices] = sep2d[keep_mask]
    match = d2d < 1
else:
    idx,d2d,_ = coo_g.match_to_catalog_sky(coo)
    match = d2d.arcsec < 1
    
g = g[match]
cols['dK'][idx[match]] = g['K_ext']-g['Kmag']
    
t1 = time.time()
print("matched to 2MASS, XSC after %.1f seconds"%(t1-t0))
t0 = time.time()

### MATCH TO LS9
cnames = {'rKmag': 'MAG_R', 'gKmag': 'MAG_G', 'zKmag': 'MAG_Z',
          'W1mag': 'MAG_W1', 'W2mag': 'MAG_W2'}
ls = "ls9_magnitudes.fits"
with fits.open(ls, memmap=True) as hdul:
    g = hdul[1].data

coo_g = SkyCoord(g['RA'],g['DEC'],frame='icrs',unit='deg')
i180 = g['RA']<180

idx1,d2d1,_ = coo_g[i180].match_to_catalog_sky(coo)
idx2,d2d2,_ = coo_g[~i180].match_to_catalog_sky(coo)
idx = np.full(len(g),0,'int64')
d2d = np.full(len(g),0,'float32')
d2d[i180] = d2d1.arcsec 
d2d[~i180] = d2d2.arcsec
idx[i180] = idx1
idx[~i180] = idx2

match = d2d < 3
g = g[match]
for c in cnames.keys():
    cname = cnames[c]
    cols[c][idx[match]] = g[cname]
    
cols['griz_ref'][idx[match]] = 3
cols['W_ref'][idx[match]] = 3

t1 = time.time()
print("matched to LS9 after %.1f seconds"%(t1-t0))


### MATCH TO DELVE
cnames = {'rKmag': 'mag_auto_r', 'gKmag': 'mag_auto_g',
          'iKmag': 'mag_auto_i', 'zKmag': 'mag_auto_z'}
de = "delve_magnitudes.fits"
with fits.open(de, memmap=True) as hdul:
    g = hdul[1].data

coo_g = SkyCoord(g['ra'],g['dec'],frame='icrs',unit='deg')
if smallcat:
    ind, sep2d, _ = match_coordinates_sky(coo, coo_g,storekdtree=False)
    sep2d = sep2d.arcsec
    # Step 2: Sort by distance (closest first)
    isort = np.argsort(sep2d)
    ind_sorted = ind[isort]
    # Step 3: Find unique a matches (closest b for each a)
    _, iu = np.unique(ind_sorted, return_index=True)
    keep_mask = np.zeros(len(ind_sorted), dtype=bool)
    keep_mask[iu] = True
    # Step 4: Reverse the sorting to original b order
    rsort = np.argsort(isort)
    keep_mask = keep_mask[rsort]
    # Step 5: Build final index array (length = len(c_a))
    idx = np.full(len(coo_g), -1, dtype=int)
    d2d = np.full(len(coo_g), np.inf)
    # For all valid b[i] with unique closest a, assign
    a_indices = ind[keep_mask]
    b_indices = np.nonzero(keep_mask)[0]
    idx[a_indices] = b_indices
    d2d[a_indices] = sep2d[keep_mask]
    match = np.logical_and(d2d < 3,
                           cols['griz_ref'][idx]!=3)
else:
    idx,d2d,_ = coo_g.match_to_catalog_sky(coo)
    match = np.logical_and(d2d.arcsec < 3,
                       cols['griz_ref'][idx]!=3)
g = g[match]
for c in cnames.keys():
    cname = cnames[c]
    cols[c][idx[match]] = g[cname]
    
cols['griz_ref'][idx[match]] = 2

t1 = time.time()
print("matched to DELVE after %.1f seconds"%(t1-t0))
t0 = time.time()



### MATCH TO LS10-South
cnames = {'rKmag': 'r_mag', 'gKmag': 'g_mag', 'iKmag': 'i_mag', 'zKmag': 'z_mag',
          'W1mag': 'W1_mag', 'W2mag': 'W2_mag'}
ls = "stacked_legacy_magnitudes.fits"
with fits.open(ls, memmap=True) as hdul:
    g = hdul[1].data

coo_g = SkyCoord(g['RA'],g['Dec'],frame='icrs',unit='deg')
if smallcat:
    ind, sep2d, _ = match_coordinates_sky(coo, coo_g, storekdtree=False)
    sep2d = sep2d.arcsec
    # Step 2: Sort by distance (closest first)
    isort = np.argsort(sep2d)
    ind_sorted = ind[isort]
    # Step 3: Find unique a matches (closest b for each a)
    _, iu = np.unique(ind_sorted, return_index=True)
    keep_mask = np.zeros(len(ind_sorted), dtype=bool)
    keep_mask[iu] = True
    # Step 4: Reverse the sorting to original b order
    rsort = np.argsort(isort)
    keep_mask = keep_mask[rsort]
    # Step 5: Build final index array (length = len(c_a))
    idx = np.full(len(coo_g), -1, dtype=int)
    d2d = np.full(len(coo_g), np.inf)
    # For all valid b[i] with unique closest a, assign
    a_indices = ind[keep_mask]
    b_indices = np.nonzero(keep_mask)[0]
    idx[a_indices] = b_indices
    d2d[a_indices] = sep2d[keep_mask]
    match = np.logical_and(d2d < 3,
                           cols['griz_ref'][idx]!=3)
else:
    idx,d2d,_ = coo_g.match_to_catalog_sky(coo)
    match = np.logical_and(d2d.arcsec < 3,
                       cols['griz_ref'][idx]!=3)
                       
g = g[match]
for c in cnames.keys():
    cname = cnames[c]
    cols[c][idx[match]] = g[cname]
    
cols['griz_ref'][idx[match]] = 4
cols['W_ref'][idx[match]] = 4

t1 = time.time()
print("matched to LS10-South after %.1f seconds"%(t1-t0))
t0 = time.time()

### MATCH TO SGA-HECATE-GLADE1
big = "../concat_hecate_sga_glade1_magnitudes.fits"
with fits.open(big, memmap=True) as hdul:
    bigg = hdul[1].data
inbig = np.any(((t['id_cat'] & (1 << 0)) != 0,
                (t['id_cat'] & (1 << 1)) != 0,
                (t['id_cat'] & (1 << 2)) != 0),
               axis=0)
cnames = {'rKmag': 'R', 'gKmag': 'G', 'iKmag': 'I', 'zKmag': 'Z'}

# SGA-HECATE
g = bigg[bigg['extmag_ref']<3] 
coo_g = SkyCoord(g['RA'],g['DEC'],frame='icrs',unit='deg')
idx,d2d,_ = coo_g.match_to_catalog_sky(coo)

match = np.logical_and(d2d.arcsec < 10, 
                       inbig[idx])

g = g[match]
for c in cnames.keys():
    cname = cnames[c]
    cols[c][idx[match]] = g[cname]

cols['griz_ref'][idx[match]] = 4 + g['extmag_ref']

# GLADE1
g = bigg[bigg['extmag_ref']==3] 
coo_g = SkyCoord(g['RA'],g['DEC'],frame='icrs',unit='deg')
idx,d2d,_ = coo_g.match_to_catalog_sky(coo)

match = np.all((d2d.arcsec < 10,
                cols['griz_ref'][idx]<2,
                inbig[idx],
                np.logical_or(np.isnan(cols['gKmag'][idx]),
                              cols['gKmag'][idx]>g['G']+2),
                ),axis=0)

g = g[match]
bgk = g['G'] - cols['gKmag'][idx[match]]

cols['gKmag'][idx[match]] = g['G']

for c in ['rKmag','iKmag','zKmag']:
    cols[c][idx[match]] += bgk
    
cols['griz_ref'][idx[match]] = 4 + g['extmag_ref']

t1 = time.time()
print("matched to big catalogs after %.1f seconds"%(t1-t0))
t0 = time.time()

Cols = []
for name in t.columns.names:
    print(name)
    data = t[name]
    dtype = data.dtype
    
    if np.issubdtype(dtype, np.floating):
        format_code = 'E'  # 32-bit float
    elif dtype==np.dtype('int8'):
        format_code = 'B'  # 8-bit integer
    elif dtype==np.dtype('int16'):
        format_code = 'I'  # 16-bit integer
    elif np.issubdtype(dtype, np.integer):
        format_code = 'J'  # 32-bit integer
    else:
        raise TypeError(f"Unsupported column type for '{name}': {dtype}")

    col = fits.Column(name=name, array=data, format=format_code)
    Cols.append(col)
    
for name in list(cols.keys()):
    data = cols[name]
    print(name)
    dtype = data.dtype
    
    if np.issubdtype(dtype, np.floating):
        format_code = 'E'  # 32-bit float
    elif dtype==np.dtype('int8'):
        format_code = 'B'  # 8-bit integer
    elif dtype==np.dtype('int16'):
        format_code = 'I'  # 16-bit integer
    elif np.issubdtype(dtype, np.integer):
        format_code = 'J'  # 32-bit integer
    else:
        raise TypeError(f"Unsupported column type for '{name}': {dtype}")

    col = fits.Column(name=name, array=data, format=format_code)
    Cols.append(col)
    
hdu = fits.BinTableHDU.from_columns(Cols)


print("defined new hdu","(%d)"%(time.time()-t0))
hdu.writeto(cat_out, overwrite=True)
print("finished task","(%d)"%(time.time()-t0))