# SCRR
This repository contains code and scripts to accompany an article in Biostatistics "Shifting-corrected regularized regression for 1H NMR metabolomics identification and quantification". 

The file source_functions.R contains all user-defined functions supporting the method proposed in the paper.

The folder Input includes preprocessed mixture data and three reference libraries discussed in Section 5 of the paper. Specifically, ExperimentalMixture.RData contains preprocessed NMR spectra of the three experimental mixtures described in Section 5.1. Serum.RData consists of three preprocessed NMR replications of breast cancer serum samples used in Section 5.2. Finally, Reference61.RData, Reference101.RData, Reference104.RData, and Reference200.RData are spectra libraries of size 61, 101, 104, and 200, respectively, as mentioned in Section 5 of the paper.

The scripts Simulation1.R and Simulation2.R can be used to reproduce Simulation 1 and Simulation 2 in Section 4 of the paper, respectively.

The script RealData.R can be used to reproduce the real data analysis described in Section 5 of the paper.
