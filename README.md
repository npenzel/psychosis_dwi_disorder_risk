Psychosis DWI Disorder Risk
===========================

This repository contains code and documentation for the analyses presented 
in our study on diffusion-weighted imaging (DWI) and psychosis risk, using 
simulated data and canonical correlation analysis (CCA).


How to Use
----------

Before running any scripts, set your working directory to the folder where 
you've saved this repository, e.g.:

    setwd("path/to/psychosis_dwi_disorder_risk")


Main Result Files
-----------------

- pronia_psychosis_risk_dwi_paper.html
- pronia_psychosis_cvresults_paper.html

  These two HTML files show the final results used in the publication.


Analysis Scripts (Rmd and HTML)
-------------------------------

- pronia_psychosis_risk_dwi_simulation.Rmd / .html
- pronia_psychosis_cvresults_simulation.Rmd / .html

  These contain the full pipeline to reproduce the analyses step by step.


Cross-Validation
----------------

- pronia_psychosis_risk_dwi_crossvalidate_simulation.R

  This script runs the entire cross-validation pipeline on the simulated data.


Core Functions
--------------

The following R scripts define reusable functions for CCA, bootstrapping, 
permutation testing, and preprocessing:

- bootstrap_cca_function.R
- cca_function.R
- center_function.R
- permute_cca_function.R
- permute_dataset.R


Folder Overview
---------------

cv_function/
    Contains all scripts related to CCA cross-validation as used in 
    pronia_psychosis_risk_dwi_crossvalidate_simulation.R

simulated_data/
    Includes:
      - Simulated datasets
      - Permutation indices (from permute_dataset.R)
      - Intermediate result files
      - simulate_date.R: script to generate new simulated data


