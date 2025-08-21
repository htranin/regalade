# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 11:58:12 2025

@author: hgtra
"""

import numpy as np
from astropy.table import Table
from astropy.io import fits
import os
from tqdm import tqdm
import time
t0=time.time()
home = os.path.expanduser("~")
os.chdir(home+"/Downloads/Catalogs/regalade")

# Parameters
input_path = "master_catalog_thin.fits"
output_path = "master_catalog_thin_without_duplicates.fits"
ra_chunk_size = 10  # degrees
tmatch_cols = [f"t{i}_match_idx" for i in range(16)]

def priority_score(row):
    """Higher is better: more ndist, brighter mags."""
    score = row['n_dist'] * 1e6  + row['R1']# prioritize ndist, R1
    return score

def process_chunk(table_chunk):
    # Keep track of rows that are modified
    modified_mask = np.zeros(len(table_chunk), dtype=bool)

    for tcol in tmatch_cols:
        match_ids = table_chunk[tcol]
        valid = match_ids > 0
        if not np.any(valid):
            continue

        # Map match_id -> row indices
        match_map = {}
        for idx, mid in enumerate(match_ids):
            if mid <= 0:
                continue
            match_map.setdefault(mid, []).append(idx)

        # Process duplicates
        for mid, idx_list in match_map.items():
            if len(idx_list) < 2:
                continue

            # Compute priority scores
            scores = [priority_score(table_chunk[i]) for i in idx_list]
            best_idx = idx_list[np.argmax(scores)]

            for i in idx_list:
                if i != best_idx:
                    table_chunk[tcol][i] = 0  # Remove association
                    modified_mask[i] = True

    # Return only rows that were changed
    return table_chunk[modified_mask]

# Process the table in RA chunks
with fits.open(input_path, memmap=True) as hdul:
    full_table = Table(hdul[1].data)
    print("read table after %d seconds"%(time.time()-t0))
    score = full_table['n_dist'] * 1e6  + full_table['R1']
    full_table = full_table[np.argsort(score)[::-1]]
    print("sorted table after %d seconds"%(time.time()-t0))
    all_chunks = []
    l = len(full_table)
    for ra_min in tqdm(range(0, 360, ra_chunk_size)):
        ra_max = ra_min + ra_chunk_size        
        mask = (full_table['ra'] >= ra_min) & (full_table['ra'] < ra_max)
        chunk = full_table[mask]
        if len(chunk) == 0:
            continue
        chunk['id'] = np.arange(len(chunk))
        chunk.write("test.fits",overwrite=True)
        print(f"Processing RA chunk {ra_min}–{ra_max} with {len(chunk)} rows...")
        cmd = (
            'java -jar C:/Users/hgtra/OneDrive/Documents/stilts.jar '
            'tmatch1 matcher="sky" params="3" in="test.fits" '
            'values="ra dec" action="keep1" out="test_clean.fits"'
        ) 
        os.system(cmd)
        cleaned_chunk = Table.read("test_clean.fits")
        no_duplicate_chunk = chunk[np.isin(chunk['id'],cleaned_chunk['id'])]
        print("identified %d no-duplicates"%len(no_duplicate_chunk))
        all_chunks.append(no_duplicate_chunk)

# Combine and save
print("Stacking all chunks...")
final_table = Table(np.hstack([chunk.as_array() for chunk in all_chunks]))
final_table.write(output_path, overwrite=True)
print(f"Saved cleaned table to {output_path}")
