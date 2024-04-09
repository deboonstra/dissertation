# meanstr

This sub-directory of `scripts` contains script files focused on exploring the mean structure selection properties of different model selection methods. The descriptions of the script files are below.

`norm_aic_sim.R`

compares the performance of AIC and QIC for a *normally* distributed response vector. This comparison can only occur if the working correlation structure is independence. The results of this simulation study can be found at 

> `./outputs/mean/norm-aic-sim`,

where the specific results are stored as `norm_aic_sim.RData`.

`norm_gls_bqicu_sim.R`

runs the same simulation as `norm-gls-sim`; however, the GOF term in QICu and BQICu were adjusted by only transforming the *normally* distributed response vector $\mathbf{y}$ with $\mathbf{V}_{r}$ based on the unstructured correlation structure. One should see the `2023-05-24` notes for more detail. The results of this simulation study can be found at 

> `./outputs/mean/norm-gls-bqicu-sim`,

where the specific results are stored as `norm_gls_bqicu_sim.RData`.

`norm_gls_sim.R`

runs a simulation to determine the effectiveness of using the GLS approach with GEEs to better determine the GOF term in QIC by only transforming the *normally* distributed response vector $\mathbf{y}$ with $\mathbf{V}_{r}$ based on the unstructured correlation structure. One should see the `2023-04-20` and `2023-05-18` notes for more detail. The results of this simulation study can be found at 

> `./outputs/mean/norm-gls-sim`,

where the specific results are stored as `norm_gls_sim.RData`.

`norm_glsx_aic_sim.R`

runs a simulation to determine the effectiveness of using the GLS approach with GEE to transform the *normally* distributed response vector ($\mathbf{y}$) and covariate matrix ($\mathbf{X}$) by $\mathbf{V}_{r}$ based on the unstructured correlation structure. Then, uses AIC and BIC to select the mean structure. One should see the `2023-06-08` notes for more information. The results of this simulation study can be found at 

> `./outputs/mean/norm-glsx-aic-sim`,

where the specific results are stored as `norm_glsx_aic_sim.RData`.

`norm_glsx_beta_sim.R`

examines the relationship between the mean parameter estimates based on fitting GLMs via GEEs and GLMs after transforming the *normally* distributed response vector ($\mathbf{y}$) and covariate matrix ($\mathbf{X}$) by $\mathbf{V}_{r}$ based on the unstructured correlation structure. The results of this simulation study can be found at 

> `./outputs/mean/norm-glsx-beta-sim`,

where the specific results are stored as `norm_glsx_beta_sim.RData`.

`norm_glsx_bqicu_sim.R`

runs a simulation to determine the effectiveness of using the GLS approach with GEE, where the GOF term in QIC, QICu and BQICu were adjusted by transforming the *normally* distributed response vector ($\mathbf{y}$) and covariate matrix ($\mathbf{X}$) by $\mathbf{V}_{r}$ based on the unstructured correlation structure. One should see the `2023-05-31` notes for more information. The results of this simulation study can be found at 

> `./outputs/mean/norm-glsx-bqicu-sim`,

where the specific results are stored as `norm_glsx_bqicu_sim.RData`.