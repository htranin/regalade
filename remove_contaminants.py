# -*- coding: utf-8 -*-
"""
Created on Wed Apr 23 09:56:21 2025

@author: hgtra
"""


from astropy.io import fits
import os
from astropy.coordinates import SkyCoord
import numpy as np
from astropy.table import Table
from tqdm import tqdm
import time

t0= time.time()

# takes about 25 minutes to run
home = os.path.expanduser("~")
os.chdir(home+"/Downloads/Catalogs/regalade")

cat1 = "master_catalog_thin_enhanced.fits"



os.chdir("../../../Downloads/Catalogs/regalade")

cat_out2 = cat1[:-5]+'_without_contam.fits'
cat_out = cat1[:-5]+'_contam.fits'
print(cat_out)


tmatch_cols = [f"t{i}_match_idx" for i in range(1, 17)][:]


with fits.open(cat1, memmap=True) as hdul:
    
    t1 = hdul[1].data  # Access the table
    
    print("catalog read","(%d)"%(time.time()-t0))
    
    idcat = t1['id_cat'].astype(np.uint32)
    ncat = np.unpackbits(idcat[:, None].view(np.uint8), axis=1).sum(axis=1)
    print("defined ncat",ncat.shape,"(%d)"%(time.time()-t0))
    
    if 1:
        mask_orphans = ncat<=2
        coo = SkyCoord(t1['ra'],t1['dec'],frame='icrs',unit='deg')
        print("defined coordinates","(%d)"%(time.time()-t0))
        c1 = coo[mask_orphans]
        r1s = np.geomspace(20,20480,6) # steps of x4
        mask = np.full(len(c1),0,dtype='bool')
        
        for i in tqdm(range(len(r1s)-1)):
            r1min = r1s[i]
            r1max = r1s[i+1]
            mask_r1 = np.logical_and(~mask_orphans,
                                     np.logical_and(t1['R1']>r1min,t1['R1']<=r1max))
            if np.sum(mask_r1)==0:
                continue
            c2 = coo[mask_r1]
            idx,d2d,_ = c1.match_to_catalog_sky(c2)
            
            # Compute position angle difference
            posang_b = c1.position_angle(c2[idx]).deg - np.asarray(t1[mask_r1][idx]['PA'])
    
            # Apply elliptical containment condition
            R1_b = np.asarray(t1[mask_r1][idx]['R1'])
            R2_b = np.asarray(t1[mask_r1][idx]['R2'])
            
            factor = 1
            rad_min = factor * (R1_b * R2_b / np.sqrt(R1_b**2 * np.sin(np.radians(posang_b))**2 + R2_b**2 * np.cos(np.radians(posang_b))**2))
            match = d2d.arcsec < rad_min
            print("matched",np.sum(match),"objects on a total of",len(match))
            mask = np.logical_or(mask,match)
    
        # add MW internal matches
        mw_orphans = np.logical_and(ncat==1,np.abs(coo.galactic.b.deg)<np.abs(15*np.cos(coo.galactic.l.deg/360*np.pi)**2)+1)
        
        print("original orphans:",np.sum(mask_orphans),"(%d)"%(time.time()-t0))
        mask_orphans[np.where(mask_orphans)[0][~mask]] = False
        print("removed orphans:",np.sum(mask_orphans),"(%d)"%(time.time()-t0))
        
    
    
    
    # NOW, FILTER OUT CONTAMINANTS
    mask_orphans = np.logical_or(mask_orphans,np.isnan(t1['D_tmean']))
    
    t13idx =  ((t1['id_cat'] & (1 << 9)) != 0)
    rmag = t1['rKmag']
    r1 = t1['R1']
    sdss_orp = (t1['id_cat'] == 2**12)
    not_primary = np.abs(t1['ref_D_input']-5.5)>2
    
        
    bogus_contam = np.all((t1['W_ref']==0,
                           t1['griz_ref']==0,
                           ncat==1,
                           t1['id_cat'] != 2**2, #HECATE
                           t1['id_cat'] != 2**5, #Cosmicflows                               
                           ),axis=0)
    
    desipv = (t1['id_cat'] == 2**3)
    td = Table(t1[desipv])
    td = td[np.argsort(td['R1'])[::-1]] # sort by decreasing R1
    td.write("test.fits",overwrite=True)
    cmd = (
        'java -jar C:/Users/hgtra/OneDrive/Documents/stilts.jar '
        'tmatch1 matcher="sky" params="120" in="test.fits" '
        'values="ra dec" action="keep1" out="test_clean.fits"'
    ) # keep 1 per group. If bogus, will be removed anyway
    os.system(cmd)
    td2 = Table.read("test_clean.fits")
    co1 = SkyCoord(td['ra'],td['dec'],frame='icrs',unit='deg')
    co2 = SkyCoord(td2['ra'],td2['dec'],frame='icrs',unit='deg')
    _,d2d,_ = co1.match_to_catalog_sky(co2)
    tmatch = d2d.arcsec<1
    desipv[np.where(desipv)[0][tmatch]] = False # we keep ungrouped orphans


    bogus_contam = np.logical_or(bogus_contam,
                                 desipv)
    
    
    print("identified %d boguses"%np.sum(bogus_contam),"(%d)"%(time.time()-t0))
    

    
    cut = np.full(len(t1),0.175)
    gbp = t1['G']-t1['BP']
    g = t1['G']
    logsep = np.log10(t1['angDist'])
    
    cut[t13idx] += 0.225
    # remove stars in Gaia angDist-color plane
    gaia_stars = np.all((gbp < np.minimum(logsep+1,-2*logsep,-10001*logsep+10000*cut),
                         ncat<4,
                         not_primary, # primary dist
                         ),axis=0)
    delve_ls = np.abs(t1['griz_ref']-3)<=1
    pa = t1['pa']
    
    # remove stars in R1-rK plane
    stars = np.logical_or(gaia_stars,
                          r1<(ncat<4)*10**((24-rmag)/10)*delve_ls*(pa!=0))
    
    # remove PM stars
    stars = np.logical_or(stars,
                          np.all(((logsep<0)*(gbp<0)*t1['PM']>1,
                                 (t1['id_cat']  & (1 << 2)) == 0, #HECATE
                                 not_primary,
                                 ),axis=0))
    
    desipv = np.logical_or((t1['id_cat'] & (1 << 3)) != 0,
                           (t1['id_cat'] & (1 << 4)) != 0,)
    print(np.sum(desipv),"DESI matches")
    
    print("identified %d stars"%np.sum(stars),"(%d)"%(time.time()-t0))
    
    ## TODO: ADD MATCH TO DUNCAN+2022 HERE. + ~20-25 min CDS xmatch
    # notcontam = np.nanmax([match['pstar'], match['pstar_cds']], axis=0) < 0.5
    
    mask = np.any((mask_orphans,
                   mw_orphans,
                   stars,
                   bogus_contam),axis=0)
    
    t1['contam_type'][mask] = 0
    t1['contam_type'][mask_orphans] += 1
    t1['contam_type'][stars] += 2
    t1['contam_type'][bogus_contam] += 4
    t1['contam_type'][mw_orphans] += 8
    
    print("replaced contam_type column","(%d)"%(time.time()-t0))
    filtered_data = t1[mask]
    
    hdul[1] = fits.BinTableHDU(data=filtered_data, header=hdul[1].header)
    print("defined new hdu","(%d)"%(time.time()-t0))
    
    hdul.writeto(cat_out, overwrite=True)
    print("finished task","(%d)"%(time.time()-t0))
    
    filtered_data = t1[~mask]

    hdul[1] = fits.BinTableHDU(data=filtered_data, header=hdul[1].header)
    print("defined new hdu","(%d)"%(time.time()-t0))
    hdul.writeto(cat_out2, overwrite=True)
    print("finished task","(%d)"%(time.time()-t0))
    
    