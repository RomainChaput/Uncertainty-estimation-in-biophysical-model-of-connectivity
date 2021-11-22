# Codes for Uncertainty propagation into numerical biophysical model of fish larval connectivity in the Florida Keys

This repository contains the post-processing codes necessary for the uncertainty analysis performed in Chaput R., Sochala P., Miron P., Kourafalou V.H., Iskandarani M., 2022. Uncertainty propagation into numerical biophysical model of fish larval connectivity in the Florida Keys. Submitted to ICES Journal of Marine Science.

In this study, we investigate the impacts of seven uncertain biological parameters on simulated larval connectivity in the Florida Keys using polynomial chaos surrogate. These parameters describe biological traits and behaviors --such as mortality, swimming abilities and orientation-- and modulate larval settlement as well as dispersal forecasts. However, these parameters are poorly constrained by observations and vary naturally between individual larvae. In the study, we characterize the input uncertainties with probability density functions informed by previous studies of Abudefduf saxatilis. The parametric domain is sampled via ensemble calculations, then a polynomial-based surrogate is built to explicitly approximate the dependence of the model outputs on the uncertain model inputs which enables a robust statistical analysis of uncertainties. This approach allows the computation of probabilistic dispersal kernels that are further analyzed to understand the impact of the parameter uncertainties. 

We share here a Matlab code that allows to replicate the method presented in the publication: computation of the polynomial chaos coefficients, construction of the model surrogate, estimation of the dispersal kernels, validation of the surrogate, and Sobol indices. This code allows to plot some of the figures presented in the publication for the reef of the Lower Keys. Unfortunately, we cannot share all of the CMS outputs due to contraints in size. This code, therefore, starts with the dispersal kernels pre-computed for the CMS samples(352 dispersal kernels for the Lower Keys, one for each combination of input parameters values) and the 100 validation runs (100 dispersal kernels for the Lower Keys, one for each validation combination) 
The surrogates computed here are specific to our study: we considered 5 uncertain input parameters, uniformly distributed on their ranges, and used a sparse matrix to identify our sampling points.

The biophysical model and specific orientation module used to obtain the samples are available at https://github.com/RomainChaput/connectivity-modeling-system 
