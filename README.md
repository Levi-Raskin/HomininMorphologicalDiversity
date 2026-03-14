This GitHub repository contains:
 - C++ code needed to conduct LDDMM simulations
 - R code needed to analyze the results of these simulations
 - Figures from the Raskin et al. (2026) presentation at the PaleoAnthro Society Meetings in Denver

It does not contain the simulation results we analyzed as these require far too much storage. Please reach out to [Levi Raskin](mailto:levi_raskin@berkeley.edu) if you want access to either the raw simulation data or the prossessed (calculated morphospace volume) files and I will figure out how to enable that. This code works and should be sufficient to reproduce our analysis, but be forewarned that this code is not guaranteed in any way. 

It has 3 subfolders.

data
  - contains the trees sampled from the phylogenetic posterior distribution. This is the same sample used in https://www.biorxiv.org/content/10.1101/2025.10.31.685754v1

scripts
  - Figures.R
  - AnalysisV2.R
      - Calculates volume occupied by known/simulated hominin in PCA space
  - SimulationScripts/
      - C++ files needed to perform the simulations of geometric morphometric data under LDDMM and to simulate additional tips on the hominin phylogenetic tree
    
figures
  - contains figures as svg files
