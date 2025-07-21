# From Risk to Disorder: Diffusion-Weighted Imaging in Psychosis

R-code for the analysis described in **From Risk to Disorder: White Matter Abnormalities in Recent-Onset Psychosis Compared to Clinical High-Risk –Cross-Sectional Findings from the PRONIA study** by *Nora Penzel, Kevin Kang Ik Cho*, Johanna Seitz-Holland, Anne Ruef, Giuseppe Cabras, Pedro Costa Klein, Linda A. Antonucci, Dominic Dwyer, Linda T. Betz, Shalaila Haas, Suheyla Cetin Karayumak, Fan Zhang, Grace Jacobs, Maria F Urquijo Castro, Ulrich Ettinger, Peter Falkai, Giulio Pergola, Rachel Upthegrove, Stefan Borgwardt, Paolo Brambilla, Rebekka Lencer, Eva Meisenzahl, Frauke Schultze-Lutter, Marlene Rosen, Theresa Lichtenstein, Lana Kambeitz-Ilankovic, Stephan Ruhrmann, Raimo R. K. Salokangas, Christos Pantelis, Stephen J. Wood, Alessandro Bertolino, Nikolaos Koutsouleris, *Ofer Pasternak, Joseph Kambeitz* and the PRONIA Consortium. XXX, 2025. 

Data for this project comes from the PRONIA study and is currently not publicly available. Details on the analytic approach can be found in the published article. In case of questions, please do not hesitate to contact us: [npenzel@mgh.harvard.edu](mailto:npenzel@mgh.harvard.edu)

---

## Table of Contents

1. [Overview](#overview)

This repository contains code for analyzing Diffusion Weighted Imaging (DWI) data from the PRONIA project along with a simulated dataset that might help navigate the main analysis steps.

The main objective of this project is to investigate the relationship between brain structure, psychosis risk and disorder. Specifically, the analyses aim to:

- Examine group effects on diffusion weighted imaging (DWI) variables.
- Compare DWI  measures across different clinical groups.
- Investigate associations with clinical and risk phenotypes.

2. [Installation](#installation)

Make sure to have the following R packages installed using install.packages('package_name'):
gtsummary, tidyverse, ggsci, broom, tidyr, ggbeeswarm, kableExtra, car, AICcmodavg, lubridate, emmeans, data4PCCAR, MASS, CCA, grid, rstatix, psych, corrplot, ggcorrplot, cowplot, DT

---

## Installation

To install **psychosis_dwi_disorder_risk**, follow these steps:

```bash
# Clone the repository
git clone https://github.com/npenzel/psychosis_dwi_disorder_risk

# Navigate to the project directory
cd psychosis_dwi_disorder_risk

--- 

## How to work with the scripts
make sure that you set the working directory to where you have saved **psychosis_dwi_disorder_risk** (setwd()).

