# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 17:31:05 2025

@author: hgtra
"""
from tqdm import tqdm
import numpy as np
from astropy.table import Table
from astropy.io import fits
import time
import os

home = os.path.expanduser("~")
os.chdir(home+"/Downloads/Catalogs/regalade")
t0=time.time()


fits_files = ['sga_input.fits', #0
              'glade1_input.fits', #1
              'hecate_input.fits', #2
              'desi_pv_input.fits', #3
              'desi_dr1_input.fits', #4
              'cosmicflows_input.fits', #5 (12)
              'nedlvs_d_input.fits', #6 (13)
              'nedlvs_z_input.fits', #7 (14)
              'nedlvs_r_input.fits', #8 (15)
              'gladep_input.fits', #9 (1st catalog with >1e7 rows)
              'desi_zou_input.fits', #10 (5)
              'panstarrs_input.fits', #11 (6)
              'sdss_input.fits', #12 (7)
              'gsc_input.fits', #13 (11)
              'ls10_wh24_input.fits', #14 (8)
              'delve_input.fits', #15 (10)
              ]

# order dist: 5, 6, 7, 4, 8, 0, 2, 3, 9, 1, 11, 13, 10, 14, 15, 12

cat_out = "master_catalog_thin.fits"
fname = "master_catalog.fits"

priority_order = [5, 6, 7, 4, # redshift ind, specz
                  8, 0, 2, 3, 9, 1, # galaxy cat
                  11, 13, 10, 14, # reliable big cat
                  15, 12, # less reliable
                  ]

if 1:
    with fits.open(fname,memmap=True) as hdul:
        full_table = hdul[1].data


D_best = np.full(len(full_table), np.nan,dtype='float32')
D_min = np.full(len(full_table), np.nan,dtype='float32')
D_max = np.full(len(full_table), np.nan,dtype='float32')
D_sum = np.full(len(full_table), 0., dtype='float32')
ncat =  np.full(len(full_table), 0, dtype='int8')
idcat =  np.full(len(full_table), 0, dtype='int32')
D_best_ref =  np.full(len(full_table), 0, dtype='int8')

trimmed_mean = np.full(len(full_table), np.nan, dtype='float32')

colnames = full_table.columns.names

# Loop through the catalogs in priority order
for i in tqdm(priority_order):
    match_col = f't{i}_match_idx'
    
    if match_col not in colnames:
        continue

    match_idx = full_table[match_col]

    # Only process rows where D_best is still nan and a match exists
    needs_update = (match_idx > 0) & np.isnan(D_best)
    matched = (match_idx > 0)
    if not np.any(needs_update):
        continue

    # Load corresponding catalog
    cat = Table.read(fits_files[i], memmap=True)
    #cat = cat[np.argsort(cat['ra'])]
    if np.ma.getmask(cat['D']).any():
        print("replacing masked values by NaNs")
        cat['D'].fill_value = np.nan
        cat['D'].mask = False
        
    cat['D'][cat['D']<=0] = np.nan
    
    # Get D values from matched rows
    idx = match_idx - 1
    matched_rows = match_idx[needs_update] - 1  # convert 1-based to 0-based index
    matched_D = cat['D'][matched_rows]
    arrmin = np.stack([D_min[matched], cat[idx[matched]]['D']], axis=1)
    arrmax = np.stack([D_max[matched], cat[idx[matched]]['D']], axis=1)
    arrsum = np.stack([D_sum[matched], cat[idx[matched]]['D']], axis=1)
    D_min[matched] = np.nanmin(arrmin,axis=1)
    D_max[matched] = np.nanmax(arrmax,axis=1)
    D_sum[matched] = np.nansum(arrsum,axis=1)
    ncat[matched] += (~np.isnan(cat[idx[matched]]['D'])).astype(int)
    idcat[matched] += 2**i
    # Warning: signed 32 bit integer limited to 31 catalogs
        
    D_best[needs_update] = matched_D
    D_best_ref[needs_update] = i

full_table['D_input'] = D_best
full_table['ref_D_input'] = D_best_ref

valid = ncat >= 3
trimmed_mean[valid] = (D_sum[valid] - D_min[valid] - D_max[valid]) / (ncat[valid] - 2)
trimmed_mean[~valid] = D_sum[~valid] / ncat[~valid]
full_table['D_tmean'] = trimmed_mean
full_table['n_dist'] = ncat
full_table['id_cat'] = idcat

# Apply threshold
ind = full_table['D_tmean']<2000
full_table = full_table[ind]
print(len(ind),"->",len(full_table))

cols = []
for name in full_table.columns.names:
    print(name)
    if "idx" in name :
        continue
    data = full_table[name]
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
hdu.writeto(cat_out, overwrite=True)
print("finished task","(%d)"%(time.time()-t0))

