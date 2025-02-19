# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 20:49:51 2024

@author: hgtra
"""

import numpy as np
import time
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
from astropy import units as u
import os
import subprocess
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pandas as pd
from tqdm import tqdm
from astropy.cosmology import Planck18
import astropy.cosmology.units as cu
from astropy.cosmology.units import redshift_distance
import warnings

warnings.filterwarnings('ignore')

rcParams.update({'font.size':17})

plot = 1
n = 100

# Define the distance threshold (in Mpc, example value)
distance_threshold = 2000  # Adjust as needed
ramins = np.linspace(0,360,10)

fileout = "regalade_v1.3.fits"
fileoutcsv = fileout.replace("fits","csv")


os.chdir("../../../Downloads/Catalogs/regalade")


# Define catalogs and their file paths/catalog names (placeholder paths)
catalogs = {
    "NEDD": "NEDLVS_20250128_dist.fits",
    "NEDSPZ": "NEDLVS_20250128_spec.fits",
    "GWGC": "GWGC.fits",
    "GLADE+": "GLADE+.fits",
    "GLADE1": "GLADE1.fits",
    "HECATE": "HECATE.fits",
    "NEDL": "NEDLVS_20250128_rest.fits",
    "DELVE": "delve2_morpho.fits",
    "PS1": "ps1_galaxies_with_morph.fits",
    "DESIPV": "DESI_Saulder.fits",
    "SDSS": "sdss16_morpho.fits",
    "SGA": "SGA_2020.fits",
    "WIS":"wiserep_unique.fits"
}


short_names = {name:name[0]+name[-1] for name in catalogs.keys()}

# Columns to keep
columns_to_keep = ["ra", "dec", "Bmag", "R1", "R2", "PA", "D"]

# Matching radius (in arcseconds)
match_radius = 10 * u.arcsec

# Load GLADE+ separately for matching
#catalog_glade_plus = Table.read(catalogs["GLADE+"], format="fits")
#catalog_glade_plus_coords = SkyCoord(ra=catalog_glade_plus['ra'], 
#                                     dec=catalog_glade_plus['dec'],unit='deg',frame='icrs')

# Function to process and filter a catalog


def process_catalog(path, catalog_name, apply_mask=True, ramin=None, ramax=None, coords=False):
    # Load catalog
    
    catalog = Table.read(path, format="fits")
    #catalog = catalog[np.argsort(catalog['ra'])]
    #if not('id' in catalog.colnames):
    #    catalog['id'] = np.arange(len(catalog)).astype(int)
    #    catalog.write(path,overwrite=True)
        
    
    #Table.read(path, format="fits")
    if ramin is not None:
        catalog = catalog[catalog['ra']>ramin]
    if ramax is not None:
        catalog = catalog[catalog['ra']<ramax]
    
    
    
    # Filter sources by distance threshold (assuming column 'distance')
    if 'D' in catalog.colnames:
        if not('z' in catalog.colnames) and 0:
            catalog['z'] = catalog['D']/4285.-(catalog['D']/6930)**1.8
            if "D_err" in catalog.colnames:
                catalog['z_err'] = catalog['D_err']/4285.
    else:
        catalog['D'] = len(catalog)*[np.nan]
        #catalog['D'] = catalog['z']*4285.+3200.*catalog['z']**2
        #print(catalog_name,"redshift min/max:",catalog['z'].min(),catalog['z'].max())
        if len(catalog)>0:
            catalog['D'] = (catalog['z']*cu.redshift).to(u.Mpc, cu.redshift_distance(Planck18, kind="luminosity")).value
            if "z_err" in catalog.colnames and 0:
                #catalog['D_err'] = catalog['z_err']*4285.+3200.*catalog['z_err']**2
                catalog['D_err'] = (catalog['z_err']*cu.redshift).to(u.Mpc, cu.redshift_distance(Planck18, kind="luminosity")).value
        
    
    # remove bad entries
    idx = np.logical_or(np.isnan(catalog['D']),
                        np.logical_and(catalog['D'] > 0, catalog['D'] < 26000))
    
    catalog = catalog[idx]
    
    idx = np.ones(len(catalog)).astype(bool)
    
    if apply_mask:
        #catalog['ref'] = len(catalog)*[catalog_name]
        idx = np.logical_and(catalog['D'] <= distance_threshold,
                             catalog['D'] > 0.04)
    if coords:
        return catalog[['ra','dec']], idx
        
    if catalog_name in ["DELVE","PS1","SDSS","GLADE1"]:
        catalog['R1'] /= 60.
        catalog['R2'] /= 60.
            

    return catalog, idx

# Process all catalogs and match iteratively
fileoutis = [fileoutcsv.replace(".csv","_%02d.csv"%(i/10)) for i in ramins[:-1]]

for i in tqdm(range(len(ramins)-1),position=1,leave=False):
    t0 = time.time()
    ramin = ramins[i]
    ramax = ramins[i+1]
    fileouticsv = fileoutis[i]
    compiled_catalogs = []
    coords = {}
    for catalog_name, catalog_path in tqdm(catalogs.items(),position=0,leave=True):
        
        catalog, idx = process_catalog(catalog_path, catalog_name, False, ramin, ramax, coords=True)
        coords[catalog_name] = SkyCoord(ra=np.asarray(catalog['ra']), dec=np.asarray(catalog['dec']),unit='deg',frame='icrs')
        
        
    for catalog_name, catalog_path in catalogs.items():
        if catalog_name == "GLADE+" and 0:
            continue
        
        catalog, idx = process_catalog(catalog_path, catalog_name, True, ramin, ramax)
        catalog = catalog[idx]
        #catalog.reset_index(drop=True,inplace=True)
        catalog_coords = coords[catalog_name][idx]
        print(catalog_name, len(catalog))
        
        catalog['id'].name = f'id_{catalog_name}'
        catalog['D'].name = f'D_{catalog_name}'
        
        if 'R1' in catalog.colnames:
            catalog['R1'].name = f'R1_{catalog_name}'
            catalog['R2'].name = f'R2_{catalog_name}'
            catalog['PA'].name = f'PA_{catalog_name}'
        elif len(catalog)>0: 
            catalog[f'R1_{catalog_name}'] = np.nan
            catalog[f'R2_{catalog_name}'] = np.nan
            catalog[f'PA_{catalog_name}'] = np.nan
            
        
        # Crossmatch with other catalogs to add alternative distances
        
        
        for other_name, other_path in tqdm(catalogs.items(),position=0,leave=True):
            #print(catalog_name, other_name)
            if other_name == catalog_name or len(catalog)==0:
                continue
            
            other_catalog, idx = process_catalog(other_path, other_name, False, ramin, ramax)
            other_coords = coords[other_name]
            
            if len(other_catalog)>0:
                
                #SkyCoord(ra=other_catalog['ra'], dec=other_catalog['dec'],unit='deg',frame='icrs')
        
                # Perform crossmatch
                
                idx, sep2d, _ = catalog_coords.match_to_catalog_sky(other_coords)
                match_mask = sep2d <= match_radius
                ind = np.where(match_mask)[0]
                idx_valid = idx[match_mask]
                sep2d_valid = sep2d[match_mask]
    
                # Get unique matches by keeping the closest one
                sort_order = np.argsort(sep2d_valid)  # Sort by separation (smallest first)
                ind_sorted = ind[sort_order]
                idx_sorted = idx_valid[sort_order]
                sep2d_sorted = sep2d_valid[sort_order]
                
                # Use np.unique to retain only the first (closest) match for each target
                _, unique_indices = np.unique(idx_sorted, return_index=True)
                
                ind_unique = ind_sorted[unique_indices]
                idx_unique = idx_sorted[unique_indices]                    
                    
                 
                # Add distance column for the other catalog
                catalog[f'D_{other_name}'] = np.nan  # Default to NaN
                #distances[match_mask] = other_catalog['D'][idx[match_mask]]
                catalog[f'D_{other_name}'][ind_unique] = other_catalog['D'][idx_unique]
                catalog[f'id_{other_name}'] = -1
                catalog[f'id_{other_name}'][ind_unique] = other_catalog['id'][idx_unique]
                
                catalog[f'R1_{other_name}'] = np.nan
                catalog[f'R2_{other_name}'] = np.nan
                catalog[f'PA_{other_name}'] = np.nan
                #r1s =  np.full(len(catalog), np.nan)  # Default to NaN
                #r2s =  np.full(len(catalog), np.nan)  # Default to NaN
                #pas =  np.full(len(catalog), np.nan)  # Default to NaN
                
                if 'R1' in other_catalog.colnames:
                    catalog[f'R1_{other_name}'][ind_unique] = other_catalog['R1'][idx_unique]
                    catalog[f'R2_{other_name}'][ind_unique] = other_catalog['R2'][idx_unique]
                    catalog[f'PA_{other_name}'][ind_unique] = other_catalog['PA'][idx_unique]
                
                    
                #catalog[f'R1_{other_name}'] = r1s
                #catalog[f'R2_{other_name}'] = r2s
                #catalog[f'PA_{other_name}'] = pas
            else:
                catalog[f'D_{other_name}'] = np.nan
                catalog[f'R1_{other_name}'] = np.nan
                catalog[f'R2_{other_name}'] = np.nan
                catalog[f'PA_{other_name}'] = np.nan
                catalog['id_{other_name}'] = -1
                
        if len(catalog)>0:
            
            compiled_catalogs.append(catalog.to_pandas())
            
    
    
    t1 = time.time()
    print("Elapsed time:",t1-t0,"seconds")
    # Combine all processed catalogs
    combined_catalog = pd.concat(compiled_catalogs, axis=0, ignore_index=1, join="inner")
    combined_catalog.reset_index(drop=True,inplace=True)

    
    # Adding distance columns (dist_max, dist_min, dist_median)
    dist_columns = [f'D_{name}' for name in catalogs.keys()]
    id_columns = [f'id_{name}' for name in catalogs.keys()]
    
    # Stack the column data into a 2D array for comparison
    data_stack = np.array([combined_catalog[col] for col in dist_columns])
    m = np.ma.masked_array(data_stack, np.isnan(data_stack))
    
    # Calculate the row-wise minimum while ignoring NaN values
    combined_catalog['dist_min'] = np.nanmin(m, axis=0)
    combined_catalog['dist_max'] = np.nanmax(m, axis=0)
    combined_catalog['dist_median'] = np.nanmedian(m, axis=0)
    
    
    r1_columns = [f'R1_{name}' for name in catalogs.keys()]
    data_stack = np.array([combined_catalog[col] for col in r1_columns])
    m = np.ma.masked_array(data_stack, np.isnan(data_stack))
    iR1 = np.nanargmax(m, axis=0)
    combined_catalog['R1'] = data_stack[iR1, np.arange(data_stack.shape[1])]
    r2_columns = [f'R2_{name}' for name in catalogs.keys()]
    data_stack = np.array([combined_catalog[col] for col in r2_columns])
    combined_catalog['R2'] = data_stack[iR1, np.arange(data_stack.shape[1])]#np.nanmax(data_stack, axis=0)
    pa_columns = [f'PA_{name}' for name in catalogs.keys()]
    data_stack = np.array([combined_catalog[col] for col in pa_columns])
    combined_catalog['PA'] = data_stack[iR1, np.arange(data_stack.shape[1])]#np.nanmedian(data_stack, axis=0)
    
    ref_ellipse = np.asarray(list(short_names.values()))
    ref_ellipse[0] = "  "
    combined_catalog['ref_ellipse'] = ref_ellipse[iR1]
    
    del data_stack
    
    combined_catalog.drop(r1_columns+r2_columns+pa_columns,axis=1, inplace=True)
    
    
    
    
    combined_catalog['dist_best'] = np.nan
    combined_catalog['ref_dist_best'] = "  "
    #combined_catalog['catalog_flag'] = 0
    combined_catalog['n_dist'] = 0
    for j, name in enumerate(catalogs.keys()):
        col = f'D_{name}'       
        idx2 = ~np.isnan(combined_catalog[col])            
        idx = np.logical_and(np.isnan(combined_catalog['dist_best']),
                             idx2)
        combined_catalog['dist_best'][idx] = combined_catalog[col][idx]
        combined_catalog['ref_dist_best'][idx] = short_names[name]
        #combined_catalog['catalog_flag'][idx2] += 10**j
        
        combined_catalog['n_dist'][idx2] += 1
    
    

    # convert to arcmin
    #combined_catalog['R1'][combined_catalog['ref']=='GWGC']/=60
    #combined_catalog['R2'][combined_catalog['ref']=='GWGC']/=60
    #combined_catalog['R1'][combined_catalog['ref']=='PS1']/=60
    #combined_catalog['R2'][combined_catalog['ref']=='PS1']/=60
    #combined_catalog['R1'][combined_catalog['ref']=='DELVE']/=60
    #combined_catalog['R2'][combined_catalog['ref']=='DELVE']/=60
    
    # Minimum 6 arcsec radius
    #combined_catalog['R1'][combined_catalog['R1']<0.1] = 0.1
    #combined_catalog['R2'][combined_catalog['R2']<0.1] = 0.1
    combined_catalog['R1'][np.isnan(combined_catalog['R1'])] = 0.05
    combined_catalog['R2'][np.isnan(combined_catalog['R2'])] = combined_catalog['R1'][np.isnan(combined_catalog['R2'])]
    combined_catalog['PA'][np.isnan(combined_catalog['PA'])] = 0
    
    combined_catalog.drop(combined_catalog.index[np.logical_and(combined_catalog['ref_dist_best']=='SS',combined_catalog['R1']>1./3)],inplace=True)
    combined_catalog.reset_index(drop=True,inplace=True)
    combined_catalog.drop(combined_catalog.index[np.logical_and(combined_catalog['n_dist']==2,
                                                                ~np.isnan(combined_catalog['D_GLADE+']+combined_catalog['D_GLADE1']))],inplace=True)
    
    #combined_catalog['R1'][np.ma.getmask(combined_catalog['R1'])] = 0.05
    #combined_catalog['R2'][np.ma.getmask(combined_catalog['R2'])] = combined_catalog['R1'][np.ma.getmask(combined_catalog['R2'])]
    #combined_catalog['PA'][np.ma.getmask(combined_catalog['PA'])] = 0
    #combined_catalog.to_csv(fileouticsv.replace('.csv','_a.csv'),index=False)
    combined_catalog.reset_index(drop=True,inplace=True)
    print("2",len(combined_catalog))
    df2 = combined_catalog.loc[np.minimum(combined_catalog['dist_median'],combined_catalog['dist_best'])<distance_threshold]
    #df2 = df[np.logical_and(df['ra']>ramins[k],
    #                        df['ra']<ramins[k+1])]
    df2.reset_index(drop=True,inplace=True)
    print("1",len(df2))
    df2 = df2.iloc[np.argsort(df2['n_dist'])[::-1]]
    df2.reset_index(drop=True,inplace=True)
    coo2 = SkyCoord(np.asarray(df2['ra']),np.asarray(df2['dec']),frame='icrs',unit='deg')
    if len(coo2)==0:
        continue
    
    i2, s2, _ = coo2.match_to_catalog_sky(coo2,nthneighbor=2)
    idx = s2.arcsec > 10
    df2_1 = df2.loc[idx]
    df2_1.reset_index(drop=True,inplace=True)
    coo2_1 = coo2[idx]
    df2_m = df2.loc[~idx]
    df2_m.reset_index(drop=True,inplace=True)
    coo2_m = coo2[~idx]
    original_indices = np.where(~idx)[0]
    
    for ms in tqdm(range(6)):
        i2, s2, _ = coo2_m.match_to_catalog_sky(coo2_m,nthneighbor=2)
        i3, s3, _ = coo2_m.match_to_catalog_sky(coo2_m,nthneighbor=3)
        #i4, s4, _ = coo2.match_to_catalog_sky(coo2,nthneighbor=4)
        sep_constraint_2 = s2.arcsec <= 10
        sep_constraint_3 = s3.arcsec <= 10
        #sep_constraint_4 = s4.arcsec <= ms
        keep2 = np.logical_or(np.arange(len(i2))<=i2,~np.isin(np.arange(len(i2)),i2[sep_constraint_2]))
        keep3 = np.logical_or(np.arange(len(i3))<=i3,~np.isin(np.arange(len(i3)),i3[sep_constraint_3]))
        #keep4 = np.logical_or(np.arange(len(i4))<=i4,~np.isin(np.arange(len(i4)),i4[sep_constraint_4]))
        keep = np.logical_and(keep2,keep3)
        
        original_indices = original_indices[keep]
        df2_m = df2_m.loc[keep]
        df2_m.reset_index(drop=True,inplace=True)
        coo2_m = coo2_m[keep]
    
    df2 = pd.concat([df2_1,df2_m],axis=0,ignore_index=1)
    df2.reset_index(drop=True,inplace=True)
    coo2 = coo2[np.concatenate([np.where(idx)[0],original_indices])]
        
        
    #df2 = filter_by_ra_dec_within_epsilon(df2, ra_epsilon=2e-3, dec_epsilon=2e-3)
    df2['l'] = coo2.galactic.l.deg
    df2['b'] = coo2.galactic.b.deg
    idx = np.logical_or(df2['n_dist']>1,
                        np.logical_and(df2['ref_dist_best']!='G1',
                                       np.logical_or(df2['ref_dist_best']!='P1',
                                                     np.logical_or(df2['b']>14-(np.abs(df2['l']-360*(df2['l']>180))/15)**1.5,
                                                                   df2['b']<-14+(np.abs(df2['l']-360*(df2['l']>180))/15)**1.5))))

    
    df2 = df2.loc[idx]
    coo2 = coo2[idx]
    df2.reset_index(drop=True,inplace=True)
    print(len(df2))
    df2.drop([c for c in list(df2.columns) if not(c in ['ra','dec','Bmag','R1','R2','PA','ref_ellipse','dist_min','dist_median','dist_max','dist_best','ref_dist_best','n_dist','l','b']+id_columns)],axis=1,inplace=True)
    
    
    onedist = df2['n_dist']==1
    df2_1 = df2.loc[onedist]
    df2_1.reset_index(drop=True,inplace=True)
    df2_b = df2.loc[~onedist]
    df2_b.reset_index(drop=True,inplace=True)
    coo2_1 = coo2[onedist]
    coo2_b = coo2[~onedist]
    
    sep_constraint = np.ones(1)
    r1min, r1max = df2_b['R1'].min(),df2_b['R1'].max()
    r1med = 10**((np.log10(r1min)+3*np.log10(r1max))/4)
    r1qua = 10**((np.log10(r1min)+np.log10(r1med))/2)
    
    large_r1 = df2_b['R1']>r1med
    medium_r1 = np.logical_and(~large_r1, df2_b['R1']>r1qua)
    small_r1 = df2_b['R1']<r1qua
    
    for iidx, idx in enumerate([large_r1,medium_r1,small_r1]):
        print(np.sum(idx),np.sum(np.isnan(idx)))
        if np.sum(idx)>0:
            i1,d1,_ = coo2_1.match_to_catalog_sky(coo2_b[idx])
        else:
            continue
        if np.sum(idx)>1:
            i2,d2,_ = coo2_1.match_to_catalog_sky(coo2_b[idx],nthneighbor=2)
        else:
            i2,d2 = i1,d1
        if np.sum(idx)>2:
            i3,d3,_ = coo2_1.match_to_catalog_sky(coo2_b[idx],nthneighbor=3)
        else:
            i3,d3 = i1,d1
        if np.sum(idx)>3:
            i4,d4,_ = coo2_1.match_to_catalog_sky(coo2_b[idx],nthneighbor=4)
        else:
            i4,d4 = i1,d1
        sep_constraint = np.logical_or(np.logical_or(d1.arcmin < np.asarray(df2_b[idx].iloc[i1]['R1'])*(1+int(iidx==0)),
                                                     d2.arcmin < np.asarray(df2_b[idx].iloc[i2]['R1'])*(1+int(iidx==0))),
                                       np.logical_or(d3.arcmin < np.asarray(df2_b[idx].iloc[i3]['R1'])*(1+int(iidx==0)),
                                                     d4.arcmin < np.asarray(df2_b[idx].iloc[i4]['R1'])*(1+int(iidx==0))))
        
        df2_1 = df2_1.loc[~sep_constraint]
        df2_1.reset_index(drop=True,inplace=True)
        coo2_1 = coo2_1[~sep_constraint]
        print(" ",len(df2_1))
                               
    
    df2 = pd.concat([df2_b,df2_1],axis=0,ignore_index=1)
    print(len(df2))
    
    df2['R1'][df2['R1']<=0.05] = 0.05
    df2['R2'][df2['R2']<=0.05] = 0.05    
    
    df2.reset_index(drop=True,inplace=True)
    df2 = df2.iloc[np.argsort(df2['ra'])]
    df2.reset_index(drop=True,inplace=True)
       
    df2.to_csv(fileouticsv,index=False)
    t2 = time.time()
    print("Elapsed time:",t2-t2,"seconds")
    
print("compiling ra slices...")
pd.concat([pd.read_csv(f) for f in fileoutis]).to_csv(fileoutcsv,index=False)



# Run stilts to remove duplicates
def execute(cmd):
    popen = subprocess.Popen(cmd, stdout=None, stderr=subprocess.PIPE, universal_newlines=True)
    while 1:
        stderr_char = popen.stderr.read(1)
        if stderr_char=="":
            break
        yield stderr_char 
    popen.stderr.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)



print("Compiled catalog saved.")

cmd = 'java -jar C:/Users/hgtra/OneDrive/Documents/stilts.jar tpipe in=%s ofmt=fits out=%s'%(fileoutcsv,fileout)
#cmd = 'java -jar C:/Users/hgtra/OneDrive/Documents/stilts.jar tmatch1 matcher="sky" params="10" in="%s" values="ra dec" action="keep1" out="%s" progress="time"'%(fileout,fileout_unique)
print(cmd)

for path in execute(cmd):
    print(path, end="")
    
print("Unique catalog saved.")
    
if plot:
    t = Table.read(fileout)
    #t['dist_best'] = scale_distance(t['dist_best'])
    plt.figure(figsize=(10,7))
    plt.scatter(t['dist_min'],t['dist_best'],s=1,alpha=0.02,color='k')
    x = np.array([0,1000])
    plt.plot(x,x,color='b',ls='--')
    plt.plot(x,x*3,color='r',label="y = 3x")
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.xlabel("dist_min (Mpc)")
    plt.ylabel("dist_best (Mpc)")
    plt.legend(loc=2)
    plt.savefig("regalade_distances.png")
    
    plt.figure(figsize=(10,7))
    plt.scatter(t['dist_best'][::3],t['dist_median'][::3],s=1,alpha=0.02,color='k')
    x = np.geomspace(0.01,100000,500)
    #plt.plot(x,scale_distance(x),color='y')
    plt.plot(x,x,color='b',ls='--')
    plt.plot(x,x*3.16,color='r',label="+/- 0.5 dex")
    plt.plot(x,x/3.16,color='r')
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.xlabel("dist_best (Mpc)")
    plt.ylabel("dist_median (Mpc)")
    plt.legend(loc=2)
    plt.savefig("regalade_distances_2.png")
    
    plt.figure(figsize=(8,6))
    plt.hist(t['n_dist'],rwidth=0.5,bins=range(8),align='left',color='k')
    plt.xlabel("Number of distance estimates")
    plt.ylabel("N galaxies")
    plt.gca().set_yscale('log')
    plt.tight_layout()
    plt.savefig("regalade_n_dist_log.png")
    
    plt.figure(figsize=(8,6))
    plt.hist(t['ref'],color='k',bins=range(7),rwidth=0.5,align='left')
    plt.gca().set_yscale('log')
    plt.xlabel("reference")
    # Rotating X-axis labels
    plt.xticks(rotation = 60)
    plt.ylabel("N galaxies")
    plt.tight_layout()
    plt.savefig("regalade_refs.png")
    #plt.gca().set_yscale('log')
    
    plt.figure(figsize=(8,6))
    idx = t['dist_min']>0.01
    b = np.geomspace(1,10,50)
    med = np.median(t['dist_best'][idx]/t['dist_min'][idx])    
    plt.axvline(med,label='median = %.1f'%med)
    pct90 = np.percentile(t['dist_best'][idx]/t['dist_min'][idx],90)
    plt.axvline(pct90,label='90th percentile = %.1f'%pct90,color='r')
    plt.hist(t['dist_best'][idx]/t['dist_min'][idx],color='k',bins=b)
    plt.gca().set_xscale('log')
    plt.legend(loc=1)
    plt.xlabel("dist_best / dist_min")
    plt.ylabel("N galaxies")
    plt.savefig("regalade_dist_qual.png")

    plt.figure(figsize=(8,6))
    plt.hist(t['R1'],bins=np.geomspace(0.1,1e3,100),color='k')
    plt.xlabel("R1 (arcmin)")
    plt.ylabel("N galaxies")
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.tight_layout()
    plt.savefig("regalade_R1.png")
    
