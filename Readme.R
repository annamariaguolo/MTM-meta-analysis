##  Copyright 2023 Annamaria Guolo (University of Padova) ##

Replication material for the manuscript 

Guolo A. (2023). Hierarchical multinomial processing tree models for meta-analysis of diagnostic accuracy studies. arXiv.


 1. MTM_software.R file for meta-analysis of diagnostic accuracy studies using approximate likleihood and the MTM approach;

 2. MTM_auxiliary_functions.c file needed to run the MTM estimation. Type "R CMD SHLIBMTM_auxiliary_functions.c" in your terminal to compile the file, and upload the resulting .so file in your R session, "dyn.load('MTM_auxiliary_functions.c.so')". The number of nodes and weights for the Gauss-Hermite quadrature needed for the exact likelihood and the pseudo-likelihood are fixed at 21. The number can be changed at user preference.

 3. MTM_data_application.R to replicate the data analysis included in the manuscript (Confusion assessment data from
Shi et al., Neuropsych Dis Treat, 2013).


January 2023, Padova, Italy.
