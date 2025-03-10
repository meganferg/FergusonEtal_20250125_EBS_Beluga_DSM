# FergusonEtal_20250125_EBS_Beluga_DSM
Data and code for estimating the abundance of Eastern Bering Sea (EBS) belugas using spatially explicit estimators, conventional design-based estimators, and ensemble modeling.

---
title: "Spatially Explicit Models of Eastern Bering Sea Beluga Density & Abundance"
author: "M.C. Ferguson, P.B. Conn, and J.T. Thorson"
date: "`r Sys.Date()`"
output: html_document
---

Welcome to a GitHub repository housing data and code for estimating the abundance of Eastern Bering Sea (EBS) belugas using spatially explicit estimators, conventional design-based estimators, and ensemble modeling.

This repository is associated with the paper

"Spatially explicit models of density improve estimates of Eastern Bering Sea beluga (*Delphinapterus leucas*) abundance and distribution from line-transect surveys" by M.C. Ferguson, P.B. Conn, and J.T. Thorson. (In review).

The parameter and variable estimates derived using this repository might differ slightly from the results presented in Ferguson et al. (in review), depending on which version of R and associated packages are used to run the scripts. The analysis for the manuscript was run using R v. 4.3.2.

**src**
The *FergusonEtal_20250125_EBS_Beluga_DSM/src* directory includes all .cpp code necessary to recreate the TMB models in Ferguson et al. (in review). All of the TMB models assume that the marginal likelihood of parameters, given observed beluga counts and other parameters, is a Tweedie distribution; the natural logarithmic link function is used to relate the counts to the additive predictor. The files include:
- *barrierTMB.hpp*: This file is called by spde_bnd_tw_DSM.cpp.  
- *null_tw_DSM.cpp*: Null (i.e., spatially constant) density model.
- *soap_tw_DSM.cpp*: Spatially explicit density model with a soap film smoother.
- *spde_bnd_tw_DSM.cpp*: Spatially explicit density model with an SPDE approximation to a Mat\'ern covariance function that accounts for barriers.
- *spde_tw_DSM.cpp*: Spatially explicit density model with an SPDE approximation to a Mat\'ern covariance function.
- *te_tw_DSM.cpp*: Spatially explicit density model with a tensor product smoother.
- *x_y_tw_DSM.cpp*: Spatially explicit density model with a bivariate and isotropic smoother.

**data**
The *FergusonEtal_20250125_EBS_Beluga_DSM/data* directory includes all objects necessary to recreate the models in Ferguson et al. (in review). The files include:
- *Ferguson_NSDL17_20240508_spde.Rdata*: This file contains all objects needed to run spde_tw_DSM.cpp, null_tw_DSM.cpp, and spde_tw_DSM_2017.R to build the 2017 SPDE model of EBS beluga density and estimate abundance.  
- *FergusonEtal_20250125_EBS_Beluga_DSM_data.Rdata*: This file contains all objects needed to build all spatially explicit models - except the 2017 SPDE model - of EBS beluga density and estimate abundance using the relevant scripts in the cpp and R folders.
- *FergusonEtal_20250125_EBS_Beluga_Nht_data.Rdata*: This file contains all objects needed by EBS_beluga_Nht_2022.R to estimate EBS beluga abundance in 2022 using the conventional, design-based estimator.

**inst**
The *FergusonEtal_20250125_EBS_Beluga_DSM/inst* directory includes all .R code necessary to recreate all of the EBS beluga density and abundance models in Ferguson et al. (in review). All of the spatially explicit density models assume that the marginal likelihood of parameters, given observed beluga counts and other parameters, is a Tweedie distribution; the natural logarithmic link function is used to relate the counts to the additive predictor. The files include:
- *EBS_beluga_Nht_2022.R*: Recreates the conventional, design-based estimators of EBS beluga abundance in 2022 for the full 2022 study area and the portion of the 2022 study area used to estimate abundance in 2017.
- *mcf_mod_eval_plots.R*: The function in this script is called by a number of other .R scripts to generate residual plots to help evaluate the fit of spatially explicit models to the data.
- *mgcv_spde_smooth_mcf*: These functions define the Mat\'ern SPDE model as a basis-penalty smoother using mgcv. The original script was published in Miller et al. (2020). This version was modified to identify inla.mesh class objects. (Miller, D.L., Glennie, R. & Seaton, A.E. Understanding the Stochastic Partial Differential Equation Approach to Smoothing. JABES 25, 1–16 (2020). https://doi.org/10.1007/s13253-019-00377-z).
- *NSDL1722_boot.R*: This script contains all of the code needed to run the parametric bootstrap that Ferguson et al. (in review) used to estimate uncertainty in abundance estimates from the density surface models.
- *soap_tw_DSM_2017.R*: EBS beluga abundance in 2017 using the soap film smoother, constructed in TMB and mgcv.
- *soap_tw_DSM_2022.R*: EBS beluga abundance in 2022 using the soap film smoother, constructed in TMB and mgcv.
- *spde_bnd_tw_DSM_2022.R*: EBS beluga abundance in 2022 using the SPDE approximation to a Mat\'ern covariance function that accounts for barriers. This model was constructed only in TMB.
- *spde_tw_DSM_2017.R*: EBS beluga abundance in 2017 using the SPDE approximation to a Mat\'ern covariance function, constructed in TMB and mgcv.
- *spde_tw_DSM_2022.R*: EBS beluga abundance in 2022 using the SPDE approximation to a Mat\'ern covariance function, constructed in TMB and mgcv.
- *te_tw_DSM_2017.R*: EBS beluga abundance in 2017 using the tensor product smoother, constructed in TMB and mgcv.
- *xy_tw_DSM_2017.R*: EBS beluga abundance in 2017 using the bivariate and isotropic smoother, constructed in TMB and mgcv.
