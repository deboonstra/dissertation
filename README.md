# Dissertation: Model Selection for Generalized Estimating Equations (GEEs)

## Problem areas to explore
- Create an information criteria that is better than QIC by accounting for the within-cluster correlation in the goodness-of-fit term.
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

All scripts and markdown files execute assuming the current working directory is set to the main `Dissertation` directory.

## Results
As described above all of the results from the simulations ran in `scripts` may be found in `outputs` with the appropriate file name, where each of the file names give a short description of the simulation is exploring. Below are more comprehensive descriptions of the simulations based on each file name.

- `norm-aic-sim` compares the performance of AIC and QIC. This comparison can only occur if the working correlation structure is independence.
- `norm-bqicu-sim` compares the performance of QIC, QICu, and BQICu through the traditional application of GEEs and information criterions (i.e., no transformations). One should see the `2023-05-24` notes for more detail. These simulations were updated to include selection of correlation structures and including unstructured correlation structure as an option.
- `norm-cic-sim` examines the correlation selection properties of CIC when the mean structure is properly specified via a simulation.
- `norm-gls-bqicu_sim` runs the same simulation as `norm-gls-sim`; however, the GOF term in QICu and BQICu were adjusted. One should see the `2023-05-24` notes for more detail.
- `norm-gls-sim` runs a simulation to determine the effectiveness of using the GLS approach with GEE to better determine the GOF term in QIC. One should see the `2023-04-20` and `2023-05-18` notes for more detail.
- `norm-glsx-aic-sim` runs a simulation to determine the effectiveness of using the GLS approach with GEE to transform the response and covariate matrix. Then, uses AIC and BIC to select the mean structure. One should see the `2023-06-08` notes for more information.
- `norm-glsx-beta-sim` examines the relationship between the mean parameter estimates based on fitting GLMs via GEEs and GLMs after transforming the response and design matrix by the GLS method.
- `norm-glsx-bqicu-sim` runs a simulation to determine the effectiveness of using the GLS approach with GEE, where the GOF term in QIC, QICu and BQICu were adjusted by producing a new $\mathbf{y}^{\ast}$ and $\mathbf{X}^{\ast}$. One should see the `2023-05-31` notes for more information.
- `norm-glsx-mix-qic-corstr-sim` runs a simulation to determine the effectiveness of using a mixed GLS approach with GEEs. The GOF term was a by-product of the transformed $\mathbf{y}^{\ast}$ and $\mathbf{X}^{\ast}$; however, the penalty term will be based on the original data. Unlike `norm-glsx-bqicu-sim`, joint selection of the mean structure and correlation structure were examined with these simulations. Additionally, different likelihood frameworks were used to define the GOF term.
- `norm-inv-cic-sim` investigates the correlation structure selection properties of inverted CIC, which is seen in KSD.
- `norm-inv-log-det-sim` investigates the correlation structure selection properties of the inverted logarithm of the determinants.
-  `norm-kdd-sim` investigates the correlation structure selection properties of the Kullback directed divergence (KDD). The functional form of KDD in the generalized estimating equations setting is based on the asymptotic distribution of $\hat{\boldsymbol{\beta}}.$
- `norm-ksd-sim` investigates the correlation structure selection properties of the Kullback symmetric divergence (KSD). The functional form of KSD in the generalized estimating equations setting is based on the asymptotic distribution of $\hat{\boldsymbol{\beta}}.$
- `norm-log-det-sim` investigates the correlation structure selection properties of the logarithm of the determinants, which is a matrix comparison technique that is used in the variants of the Kullback divergence.
- `norm-modified-cic-sim` examines the correlation selection properties of a modified version CIC when the mean structure is properly specified and $\boldsymbol{\Omega}_{I}$ is replaced by $\boldsymbol{\Omega}_{r}$ via a simulation. These simulations should be compared to `norm_cic_sim`.
- `norm-modified-inv-cic-sim` investigates the correlation structure selection properties of modified inverted CIC, which is seen in modified KSD.
- `norm-modified-inv-log-det-sim` investigates the correlation structure selection properties of the modified inverted logarithm of the determinants.
- `norm-modified-kdd-sim` investigates the correlation structure selection properties of the modified Kullback directed divergence (KDD~m~).
- `norm-modified-ksd-sim` investigates the correlation structure selection properties of the modified Kullback symmetric divergence (KSD~m~).
- `norm_modified-log-det-sim` investigates the correlation structure selection properties of the modified logarithm of the determinants, which is a matrix comparison technique that is used in the variants of the Kullback divergence.
- `norm-two-stage-qicu-sim` examines the mean and correlation selection properties of a two-stage model selection procedure based on CIC, QIC, and QICu.