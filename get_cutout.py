# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 17:59:38 2025

@author: hgtra
"""

from astropy.coordinates import SkyCoord
import os
import sys
import requests
from bs4 import BeautifulSoup
import shutil

os.chdir("plots")

def save_cutout(ra, dec, fname):
    """Saves a jpg cutout around ra, dec of 2.1 arcmin Field of View
    Automatically selects PanSTARRS, Legacy Survey or DSS based on coordinates
    ra: Right Ascension in degrees
    dec: Declination in degrees
    fname: filename to save the cutout
    """
    coo = SkyCoord(ra, dec, unit='deg', frame='icrs')
    b = coo.galactic.b.deg # Galactic latitude
    if dec>-30: # South is not covered by PanSTARRS
        url = f'https://ps1images.stsci.edu/cgi-bin/ps1cutouts?layer=ps1&pos={ra} {dec}&filter=color&filter=g&filter=z&filter=i&filetypes=stack&size=504&format=jpg&autoscale=99.500000&output_size=200'
        response = requests.get(url, timeout=15)
        soup = BeautifulSoup(response.text, "html.parser")
        td_tag = soup.find("td")
        img_tag = td_tag.find("img")
        url = 'https:'+img_tag['src']
    elif (abs(b)>20 # Galactic Plane were not covered by Legacy Surveys
          and (abs(dec+73)>1 or abs(ra-13)>3) # SMC center not covered
          and (abs(ra-80)>7 or abs(dec+68)>3)): # LMC center not covered
        url = 'https://www.legacysurvey.org/viewer/cutout.jpg?layer=ls-dr10&ra=%s&dec=%s&zoom=13'%(ra,dec)
    elif (abs(b)<9.5):
        url = 'https://decaps.legacysurvey.org/viewer/cutout.jpg?layer=decaps2&ra=%s&dec=%s&zoom=13'%(ra,dec)
    else: # fall back to lower resolution / depth
        url = 'https://archive.stsci.edu/cgi-bin/dss_search?v=poss2ukstu_red&r=%s&d=%s&e=J2000&h=2.1&w=2.1&f=gif'%(ra,dec)
        
        
    print(url)
    try:
        cutout = requests.get(url, stream=True, timeout=15)

        if cutout.status_code == 200:
            with open(fname, 'wb') as f:
                shutil.copyfileobj(cutout.raw, f)
            print('Image successfully Downloaded')
        else:
            print('Image Couldn\'t be retrieved')
    except requests.exceptions.Timeout:
        print('trying with dss...')
        url = 'https://archive.stsci.edu/cgi-bin/dss_search?v=poss2ukstu_red&r=%s&d=%s&e=J2000&h=2.1&w=2.1&f=gif' % (ra, dec)
        cutout = requests.get(url, stream=True, timeout=15)
    
        if cutout.status_code == 200:
            with open(fname, 'wb') as f:
                shutil.copyfileobj(cutout.raw, f)
            print('Image successfully Downloaded')
        else:
            print('Image Couldn\'t be retrieved')
            return
        
        
if __name__=="__main__":
    # optional: get command line arguments
    if len(sys.argv)>2:
        ra,dec = float(sys.argv[1]),float(sys.argv[2])
    else:
        ra,dec = 359.5929778540899, -20.3187275380

    save_cutout(ra,dec,"Be_star_cutout_example.jpg")



