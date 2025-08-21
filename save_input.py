# -*- coding: utf-8 -*-
"""
Created on Fri May  2 10:11:53 2025

@author: Hugo Tranin

Reads a fits catalog, removes duplicates 
and writes a subset of columns
"""

from astropy.io import fits
import numpy as np
import os
from astropy.coordinates import SkyCoord
from astropy.table import Table

home = os.path.expanduser("~")
os.chdir(home+"/Downloads/Catalogs/regalade")

fname = "hecate.fits"
fout = "hecate__input.fits"
keepcols = ["ra","dec","D","R1","R2","PA"]

with fits.open(fname,memmap=True) as hdul:
    
    t = hdul[1].data
    if 'R1' in t.columns.names:
        for c in keepcols:
            if np.ma.is_masked(t[c]):
                raise SystemExit("found that column %s is masked"%c)
        
        if np.nanmean(t['R1'])<3:
            print("converting arcmin to arcsec")
            t['R1']*=60
            t['R2']*=60
        if "GLADE1" in fname or "HECATE" in fname or "hyperleda" in fname:
            t['R1']*=1.3
            t['R2']*=1.3
            
        # replace nan and tiny by 3 arcsec
        idx = np.isnan(t['R1'])
        print("Fraction without R1:",np.mean(idx))
        t['R1'][np.isnan(t['R1'])] = 3.
        t['R2'][np.isnan(t['R2'])] = 3.
        t['PA'][np.isnan(t['PA'])] = 0.
        idx = t['R1']<3
        print("Fraction with tiny R1:",np.mean(idx))
        t['R2'][idx] = 3.
        t['PA'][idx] = 90.
        t['R1'][idx] = 3.
        
        
        t = t[np.argsort(np.isnan(t['D']).astype(int)*1000-t['R1'])]
   
    else:
        t = Table(t)
        t['R1'] = 3.
        t['R2'] = 3.
        t['PA'] = 0.
        for c in keepcols:
            if np.ma.is_masked(t[c]):
                raise SystemExit("found that column %s is masked"%c)
        t = t[np.argsort(t['D'])]
    
    print("read")
    t = t[np.argsort(t['ra'])]
    print("sorted")
    
    coo= SkyCoord(t['ra'],t['dec'],frame='icrs',unit='deg')
    
    # Iteratively cross-match each source with its nearest neighbor.
    # A neighbor is discarded if it lies within a factor times the 
    # elliptical radius defined by (R1, R2, PA) of either galaxy 
    N=3
    
    for i in range(N):
        idx,d2d,_=coo.match_to_catalog_sky(coo,nthneighbor=2)
        
        # Compute position angle difference
        posang_a = coo[idx].position_angle(coo).deg - np.asarray(t['PA'])
        posang_b = coo.position_angle(coo[idx]).deg - np.asarray(t[idx]['PA'])

        # Apply elliptical containment condition
        R1_a = np.asarray(t['R1'])
        R2_a = np.asarray(t['R2'])
        R1_b = np.asarray(t[idx]['R1'])
        R2_b = np.asarray(t[idx]['R2'])
        factor = 0.5
        rad_min = np.maximum(factor * np.minimum((R1_a * R2_a / np.sqrt(R1_a**2 * np.sin(np.radians(posang_a))**2 + R2_a**2 * np.cos(np.radians(posang_a))**2)),
                                                 (R1_b * R2_b / np.sqrt(R1_b**2 * np.sin(np.radians(posang_b))**2 + R2_b**2 * np.cos(np.radians(posang_b))**2))),
                             3.)
        
        match = d2d.arcsec<rad_min
        
        discard = np.logical_and(match,idx<np.arange(len(coo)))
        t = t[~discard]
        coo = coo[~discard]
        print(np.mean(discard))
        if np.mean(discard)==0:
            break
        
    t = t[np.argsort(t['ra'])]
    print("sorted")
    
    cols = []
    # Keep only RA, Dec, distance, and ellipse parameters in output.
    for name in keepcols:
        col_data = t[name].astype(np.float32)
        cols.append(fits.Column(name=name, array=col_data, format='E'))

    # Create new HDU and write
    hdul = fits.BinTableHDU.from_columns(cols)
    print("defined new hdu")
    hdul.writeto(fout, overwrite=True)
    print("finished task")

    
