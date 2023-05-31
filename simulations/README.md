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