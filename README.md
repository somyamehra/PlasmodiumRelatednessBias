This repository contains code and a minimum analysis dataset to accompany the manuscript **Systematic bias in malaria parasite relatedness estimation** by Mehra et al, 2024. This is available as a preprint on bioRxiv: <https://www.biorxiv.org/content/10.1101/2024.04.16.588675v1>.

The analysis can be reproduced by running the R Markdown files as follows:
* *Simulation_model_results.Rmd*:  this runs a simulation model of parasite ancestries, and performs numerical analyses of pairwise relatedness on simulated genotypic data.  
* *Empirical_results.Rmd*: this contains a case study of WGS data from a highly inbred parasite sample from Guyana and Colombia. 

A minimum analysis dataset, comprising a matrix of parasite genotypes (alleles encoded 0/1, missing sites encoded NA) and an accompanying data frame of variant sites, is available in the folder  *Empirical_data*.
