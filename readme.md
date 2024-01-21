# Project Title

## Overview

This repository contains the code to reproduce the analyses reported in the following publication:

**Title:** Measuring the 4 higher-order values in Schwartz’s theory: Validation of a 17-Item inventory  
**Authors:** C.M. Lechner, C. Beierlein, E. Davidov, S.H. Schwartz  
**Year:** 2024  
**Journal:** Journal of Personality Assessment  

## Steps to Reproduce Analyses

To reproduce the analyses, follow these steps:

### 1. Obtain the Data

Visit the [GESIS Panel](https://www.gesis.org/en/gesis-panel/gesis-panel-home) website and follow the steps to obtain the required data.

### 2. Download SPSS Data

Download the SPSS data distribution from the GESIS Panel and copy the files into the `01_data/spss` folder.

### 3. Generate R Data File

Run the script `00_data_preparation.R` to generate the R data file.

### 4. Reproduce Analyses
You can now reproduce all analyses by running the following scripts:

`01_analysis.R` for main analyses.
`02_analysis_mi.R` for measurement invariance analyses.


###  5. Generate Figures and Tables
Run the code chunks in preprint.Rmd for figures and tables as reported in the `preprint.Rmd` and `appendix.Rmd` for additional content reported in the supplementary online material.

### 6. Contact Information

If you encounter any issues or have questions, feel free to reach out to Clemens Lechner at [clemens.lechner@gesis.org](mailto:clemens.lechner@gesis.org).

Enjoy reproducing the analyses!
