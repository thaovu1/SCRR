# SCRR
This repository contains code and scripts to accompany an article in Biostatistics "Shifting-corrected regularized regression for 1H NMR metabolomics identification and quantification". 

The file source_functions.R contains all user-defined functions supporting the method proposed in the paper.

The scripts Simulation1.R and Simulation2.R can be used to reproduce Simulation 1 and Simulation 2 in Section 4 of the paper, respectively.

The script RealData.R can be used to reproduce the real data analysis described in Section 5 of the paper. Specifically, the folder Input includes preprocessed mixture data and four reference libraries discussed in Section 5 of the paper. Specifically, Exp_data.RData contains preprocessed NMR spectra of the three experimental mixtures described in Section 5.1. Hart_data.RData consists of three preprocessed NMR replications of breast cancer serum samples used in Section 5.2. Finally, spec_lib61.RData, spec_lib101.RData, spec_lib200.RData are spectra libraries of size 61, 101, and 200 used to evaluate experimental mixtures in Section 5.1. Furthermore, spec_lib104.RData is the reference library of size 104 used for serum samples in Hart et al. in Section 5.2.




