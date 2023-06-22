# Simulations for dissertation

## File directories
Within the `simulations` directory the following directories are orgranized as follows.

- `explore` contains any exploratory script files that may or may not execute.
- `outputs` contains all the results generated from `scripts` and will have file names associated with the script files that generated the results.
- `R` contains all `R` script files used by other scripts. To load all the scripts:
```r
R <- list.files(path = "./simulations/R", pattern = "*.R", full.names = TRUE)
sapply(R, source, .GlobalEnv)
```
- `scripts` contains all the script files used to execute certain tasks.

All scripts and markdown files execute assuming the current working directory is set to the main `Dissertation` directory.

## Results
As described above all of the results from the simulations ran in `scripts` may be found in `outputs` with the appropriate file name, where each of the file names give a short description of the simulation is exploring. Below are more comprehensive descriptions of the simulations based on each file name.

- `norm_aic_sim` compares the performance of AIC and QIC. This comparison can only occur if the working correlation structure is independence.
- `norm_bqicu_sim` compares the performance of QIC, QICu, and BQICu through the traditional application of GEEs and information criterions (i.e., no transformations). One should see the `2023_05_2023` notes for more detail. These simulations were updated to include selection of correlation structures and including unstructured correlation structure as an option.
- `norm_gls_bqicu_sim` runs the same simulation as `norm_gls_sim`; however, the GOF term in QICu and BQICu were adjusted. One should see the `2023_05_2023` notes for more detail.
- `norm_gls_sim` runs a simulation to determine the effectiveness of using the GLS approach with GEE to better determine the GOF term in QIC. One should see the `2023_04_20` and `2023_05_18` notes for more detail.
- `norm_glsx_aic_sim` runs a simulation to determine the effectiveness of using the GLS approach with GEE to transform the response and covariate matrix. Then, uses AIC and BIC to select the mean structure. One should see the `2023_06_08` notes for more information.
- `norm_glsx_bqicu_sim` runs a simulation to determine the effectiveness of using the GLS approach with GEE, where the GOF term in QIC, QICu and BQICu were adjusted by producing a new $\mathbf{y}^{\ast}$ and $\mathbf{X}^{\ast}$. One should see the `2023_05_31` notes for more information.
- `norm_glsx_mix_qic_corstr_sim` runs a simulation to determine the effectiveness of using a mixed GLS approach with GEEs. The GOF term was a by-product of the transformed $\mathbf{y}^{\ast}$ and $\mathbf{X}^{\ast}$; however, the penalty term will be based on the original data. Unlike `norm_glsx_bqicu_sim`, joint selection of the mean structure and correlation structure were examined with these simulations. Additionally, different likelihood frameworks were used to define the GOF term.