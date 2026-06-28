# REGALADE: Revised Galaxy List for the Advanced Detector Era

We present the **REGALADE catalog**, a compilation of galaxies with their distances, magnitudes, angular extents, and stellar masses, designed to support transient and multi-messenger astronomy.  

By merging data from curated galaxy catalogs and deep surveys, REGALADE provides a more complete census of galaxies in the Local Universe, enabling accurate distance estimates and improved matching for transients and gravitational-wave events.

---

## Contents of the Catalog

- Nearly **80 million galaxies** within **D < 2000 Mpc**
- Derived from widely used catalogs and surveys:
  - GLADE  
  - NED Local Volume Survey  
  - Siena Galaxy Atlas  
  - Legacy Surveys  
  - SDSS  
  - Pan-STARRS  
- Additional cross-matching removes duplicates and combines independent distance measurements.  

The paper describing REGALADE and published to *A&A* is available on [arXiv:2508.13267](https://arxiv.org/abs/2508.13267) (DOI [10.1051/0004-6361/202556896](https://doi.org/10.1051/0004-6361/202556896)), including case studies on:
- Optical transients  
- Gravitational-wave target selection  
- Ultraluminous X-ray sources  

---

## REGALADEv2

**REGALADEv2** revises the distance and redshift assignments of the original catalog and adds Simbad cross-matching for nearby galaxies. Key changes:

- **New distance priority** (D column): HECATE > DESI DR1 > NED-LVS-z > Cosmicflows > NED-LVS-D > NED-LVS-rest > unchanged
- **New redshift priority** (z column): DESI DR1 > NED-LVS-z > NED-LVS-rest > unchanged. Galaxies from distance-only catalogs (Cosmicflows, NED-LVS-D) that have no Simbad redshift have `z` set to NaN.
- **Distance-size outlier guard**: galaxies with `D < 100 / R1` (D in Mpc, R1 in arcsec) have their distance replaced with `D_tmean`.
- **Simbad cross-matching** (D < 200 Mpc, 5 arcsec radius): Simbad redshifts replace catalog z for distance-only sources; significant discrepancies (`|z − z_simbad| > 0.1·z + 0.001`) are flagged.
- **New columns**: `ref_z_in`, `match_offset`, `simbad_z`, `f_simbad_zdiscrepancy`

> **Input catalog versions**: prior to this processing, the HECATE and NED-LVS input catalogs were updated to **HECATEv2** and **NED-LVS v20250606** respectively.

The script `fix_distances.py` reproduces the v2 catalog from REGALADEv1.

The REGALADEv2 FITS file (~10 GB) is available at:
[Google Drive link](https://drive.google.com/file/d/1vVdR_KphBTL9E17xl4cadF1HYW_A6Zis/view?usp=sharing)

---

## Repository Contents

This repository contains code to reproduce the building of REGALADE and generate figures from the article.

Scripts:
- `save_input.py` – preprocess an input catalog for later ingestion  
- `match_input.py` – match input catalogs with each other, add cross-identification columns  
- `compute_distance.py` – add distance metrics (`D_tmean`, `D_input`, `n_dist`), remove `D > 2000 Mpc`, and save a slim master catalog  
- `clean_duplicates.py` – remove residual duplicates  
- `add_magnitudes.py` – match the master catalog to photometry from deep surveys and add magnitude columns  
- `remove_contaminants.py` – apply cleaning criteria (Section 3.5 of the paper) to remove stars, clumps, artifacts, etc.
- `fix_distances.py` – generate REGALADEv2 with revised distances, redshifts, and Simbad cross-matching

---

## Prerequisites

- Each input catalog must be downloaded locally (one FITS file per catalog).  
- Python 3.x with the following libraries:  
  - `astropy`  
  - `numpy`  
  - `tqdm`  
- [STILTS](http://www.star.bris.ac.uk/~mbt/stilts/) installed and available in your `PATH`.

---


## Data Access

The REGALADE catalog is available on VizieR at the following link:
https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/706/A284

An interactive visualization is also available at:
[https://blackpearl.blackgem.org/regalade.php](https://blackpearl.blackgem.org/regalade.php)

The full catalog can be downloaded as a single **FITS file (10 GB)**:
[https://drive.google.com/file/d/10CJa5xheifY03dOboTlhM_fvruNXxHor/view?usp=drive_link](https://drive.google.com/file/d/10CJa5xheifY03dOboTlhM_fvruNXxHor/view?usp=drive_link)

For easier handling, a **lighter version with fewer columns** is also provided in **10 CSV chunks (~500 MB each)** covering successive right ascension intervals:

* [https://blackpearl.blackgem.org/regalade_thin_0_36.csv](https://blackpearl.blackgem.org/regalade_thin_0_36.csv)
* [https://blackpearl.blackgem.org/regalade_thin_36_72.csv](https://blackpearl.blackgem.org/regalade_thin_36_72.csv)
* [https://blackpearl.blackgem.org/regalade_thin_72_108.csv](https://blackpearl.blackgem.org/regalade_thin_72_108.csv)
* [https://blackpearl.blackgem.org/regalade_thin_108_144.csv](https://blackpearl.blackgem.org/regalade_thin_108_144.csv)
* [https://blackpearl.blackgem.org/regalade_thin_144_180.csv](https://blackpearl.blackgem.org/regalade_thin_144_180.csv)
* [https://blackpearl.blackgem.org/regalade_thin_180_216.csv](https://blackpearl.blackgem.org/regalade_thin_180_216.csv)
* [https://blackpearl.blackgem.org/regalade_thin_216_252.csv](https://blackpearl.blackgem.org/regalade_thin_216_252.csv)
* [https://blackpearl.blackgem.org/regalade_thin_252_288.csv](https://blackpearl.blackgem.org/regalade_thin_252_288.csv)
* [https://blackpearl.blackgem.org/regalade_thin_288_324.csv](https://blackpearl.blackgem.org/regalade_thin_288_324.csv)
* [https://blackpearl.blackgem.org/regalade_thin_324_360.csv](https://blackpearl.blackgem.org/regalade_thin_324_360.csv)

Please feel free to contact me if you encounter any issues accessing or downloading the data.
