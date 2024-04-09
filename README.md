# Dissertation: Model Selection for Generalized Estimating Equations (GEEs)

## Problem areas to explore
- Create an information criteria that is better than QIC(R) by accounting for the within-cluster correlation in the goodness-of-fit term.
- Create an information criteria that performs better than CIC(R) by allowing more complex working correlation structures such as unstructured to be included as a valid working correlation structure.
- Data drive method to determine the variance estimator to use for inference of $\hat{\boldsymbol{\beta}}$ for GEEs.

## File directories
All directories and sub-directories are either one phrase or a combination of phrases separated by hyphens (`-`), while files within directories are either one phrase or a combination of phrases separated by underscores (`_`).

- `notes` contains files related to meetings or notes to share.
- `references` contains PDF files and notes of previous research in this area.
- `misc` contains random files that might be of use.
- `explore` contains any exploratory script files that may or may not execute.
  - This sub-directory is no longer tracked as it does not provide valuable information.
- `outputs` contains all the results generated from `scripts` and will have file names associated with the script files that generated the results.
- `scripts` contains all the script files used to execute certain tasks.
- `R` contains all `R` script files used by other scripts. To load all the scripts:

```r
R <- list.files(path = "./R", pattern = "*.R", full.names = TRUE)
sapply(R, source, .GlobalEnv)
```

All scripts and markdown files execute assuming the current working directory is set to the main `dissertation` directory.