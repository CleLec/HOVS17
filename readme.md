# Higher-Order Values Scale-17 (HOVS17) - Validation Study

## Overview

This repository contains the code to reproduce the analyses reported in the following publication:

**Title:** Measuring the 4 higher-order values in Schwartz’s theory: Validation of a 17-Item inventory  
**Authors:** C.M. Lechner, C. Beierlein, E. Davidov, S.H. Schwartz  
**Year:** 2024  
**Journal:** Journal of Personality Assessment  

The analyses results reported in the manuscript can be found in the folder `01_results` in the files `all_results.Rdata` and `mi_results.Rdata`, respectively.
These files hold results in various R data formats, such as lists, data frames / tibbles, or vectors. The tables and figures are generated from these data.
The steps outlined below will reproduce bot the `.Rdata` files and these tables and figures.

## Steps to Reproduce Analyses

To reproduce the analyses, follow these steps:

### 1. Obtain the Data

Visit the [GESIS Panel](https://www.gesis.org/en/gesis-panel/gesis-panel-home) website and follow the steps to obtain the required data.

### 2. Download SPSS Data

Download the SPSS data distribution from the GESIS Panel and copy the files into the `01_data/spss` folder.

### 3. Generate R Data File from SPSS distribution

Run the script `00_data_preparation.R` to generate the R data files from the GESIS Panel's SPSS data.

### 4. Reproduce Analyses
You can now reproduce all analyses by running the following scripts:

`01_analysis.R` for main analyses.
`02_analysis_mi.R` for measurement invariance analyses.

This will re-generate the results already stored in  `01_results/all_results.Rdata` and `01_results/mi_results.Rdata`.
Alternatively, you can skip this step and start with the results already stored there. In this case, proceed to step 5.

###  5. Generate Figures and Tables
Run the code chunks in preprint.Rmd for figures and tables as reported in the `preprint.Rmd` and `appendix.Rmd` for additional content reported in the supplementary online material.
Note that the final accepted manuscript underwent some changes (substantial textual changes, omitting the EFA results, and more) during the revision, such that `preprint.Rmd` is not identical with the final accepted manuscript.

## Contact Information

If you encounter any issues or have questions, feel free to reach out to Clemens Lechner at [clemens.lechner@gesis.org](mailto:clemens.lechner@gesis.org).

Enjoy reproducing the analyses!
