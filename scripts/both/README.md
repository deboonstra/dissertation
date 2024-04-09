# both

This sub-directory of `scripts` contains script files focused on exploring the mean structure and correlation structure selection properties of different model selection methods. The descriptions of the script files are below.

`norm_bquic_sim.R`

runs a simulation to compare the performance of QIC, QICu, and BQICu through the traditional application of GEEs and information criterions (i.e., no transformations) for a *normally* distributed response vector. One should see the `2023-05-24` notes for more detail. These simulations were updated to include selection of correlation structures and including unstructured correlation structure as an option. The results of this simulation study can be found at 

> `./outputs/both/norm-bqicu-sim`,

where the specific results are stored as `norm_bqicu_sim.RData`.

`norm_two_stage_qicu_sim.R`

runs a simulation to examine the mean and correlation selection properties of a two-stage model selection procedure based on CIC, QIC, and QICu for a *normally* distributed response vector. The results of this simulation study can be found at 

> `./outputs/both/norm-two-stage-qicu-sim`,

where the specific results are stored as `norm_two_stage_qicu_sim.RData`.

`norm_glsx_mix_qic_corstr_sim.R`

runs a simulation to determine the effectiveness of using a mixed GLS approach with GEEs. The GOF term was a by-product of transforming the *normally* distributed response vector ($\mathbf{y}$) and covariate matrix ($\mathbf{X}$) by $\mathbf{V}_{r}$ based on the unstructured correlation structure; however, the penalty term will be based on the original data. Unlike `norm_glsx_bqicu_sim.R`, joint selection of the mean structure and correlation structure were examined with these simulations. Additionally, different likelihood frameworks were used to define the GOF term.  The results of this simulation study can be found at 

> `./outputs/mean/norm-glsx-mix-qic-corstr-sim`,

where the specific results are stored as `norm_glsx_mix_qic_corstr_sim.RData`.