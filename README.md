# REGALADE: Revised Galaxy List for the Advanced Detector Era

We present the **REGALADE catalog**, a compilation of galaxies with their distances, magnitudes, angular extents, and stellar masses, designed to support transient and multi-messenger astronomy.  

By merging data from widely used galaxy catalogs and deep surveys, REGALADE provides a more complete census of galaxies in the Local Universe, enabling accurate distance estimates and improved matching for transients and gravitational-wave events.

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

A paper describing REGALADE has been submitted to *A&A* and is available on [arXiv:2508.13267](https://arxiv.org/abs/2508.13267), including case studies on:
- Optical transients  
- Gravitational-wave target selection  
- Ultraluminous X-ray sources  

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

---

## Prerequisites

- Each input catalog must be downloaded locally (one FITS file per catalog).  
- Python 3.x with the following libraries:  
  - `astropy`  
  - `numpy`  
  - `tqdm`  
- [STILTS](http://www.star.bris.ac.uk/~mbt/stilts/) installed and available in your `PATH`.  

