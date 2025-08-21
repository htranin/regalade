# Revised Galaxy List for the Advanced Detector Era

We present the REGALADE catalog, a compilation of galaxies with their distances, magnitudes, angular extent and stellar masses, designed to support transient and multi-messenger science. By merging data from a curated galaxy catalogs and wide range of deep surveys, REGALADE aims to provide a more complete census of galaxies, enabling accurate distance estimates and transient / GW event matching.

REGALADE combines data from widely used galaxy catalogs and imaging surveys, including GLADE, NED Local Volume Survey, Siena Galaxy Atlas, the Legacy Surveys, SDSS, and Pan-STARRS, with additional cross-matching to remove duplicates and combine independent distance measurements.

The REGALADE catalog contains nearly 80 million galaxies in the Local Universe (D < 2000 Mpc). A paper explaining the catalog and its contents has been submitted to A&A and is available on arxiv ([https://arxiv.org/abs/2508.13267](https://arxiv.org/abs/2508.13267)), also showcasing its use through 3 case studies: optical transients, gravitational wave target selection and ultraluminous X-ray sources.

In this repository you will find codes to reproduce the building of REGALADE and some figures of the article.

* save_input.py: preprocess an input catalog for later ingestion
* match_input.py: match input catalogs with each other, add cross-identification columns
* compute_distance.py: use cross-identification columns to add D_tmean, D_input, n_dist. Remove D > 2000 entries and save a slim master_catalog version
* clean_duplicates.py: remove residual duplicates, if any
* add_magnitudes.py: match master_catalog to photometry from deep surveys, add magnitude columns
* remove_contaminants.py: apply cleaning criteria (Section 3.5 of the article) to remove stars, galaxy clumps, artifacts etc.

Prerequisites
-------------
- Each input catalog must be downloaded locally (one FITS file per catalog)
- Python 3.x with the following libraries installed:
    * astropy
    * numpy
    * tqdm 
- STILTS installed and available in your PATH 
