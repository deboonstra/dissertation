# The function of this script file is to load the required packages to run 
# simulations.
suppressMessages(
  expr = {
    library(dplyr)
    library(Matrix)
    library(expm)
    library(geepack)
    library(gee)
    library(mvtnorm)
    library(knitr)
    library(kableExtra)
    library(rmarkdown)
    library(stats)
    library(utils)
  }
)