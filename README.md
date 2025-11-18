This directory contains the code and data needed to reproduce the results from the manuscript: _Dispersal modifies temperature dependence of ecosystem
productivity after an extreme thermal fluctuation._ We ran an outdoor freshwater mesocosm experiment in the summer of 2018 to test the hypothesis that different dispersal levels modifies the thermal sensitivity of Gross Primary Productivity in aquatic communities. 



**"Data" folder**

This folder contains data used in the manuscript. Some data files are explicitly called in the R scripts, while others contain extra supplementary data. The most important files are:

"metacom_exp_data_final.csv" : Master datasheet for autotroph and zooplankton biomass, GPP, total Nitrogen, and phytoplankton community size structure. For GLMMs/SEM. 

"phytoplankton_diversity_data.csv" : Phytoplankton community data (estimated with visual key-based ID of flow cytometry images)

"Bacterial_16S_data.csv" : Bacterial community data (estimated with 16S)

"Zoop_data.csv" : Zooplankton community data (estimated by microscopy)



**"Scripts" folder**

R scripts are roughly organized by analysis and data type and in the order presented in the manuscript. 

SEM_analyses_Fig2.R runs week-specific Structural Equation Modelling with the brms R package, and reproduces Fig 2a-d in the manuscript.
Hypothesis1_GLMMS.R runs week-specific GLMMs to test the hypothesis that the slope describing the temperature vs GPP relationship is modified by dispersal. The model is described in the Methods section of the manuscript.
Hypothesis2_Figs4-6.R generates manuscript figures4-6, which include the time series of GPP and biomass, and multivariate analyses (PERMANOVA, PERMDISP) summarizing experimental treatment effects on community size structure and taxonomic composition. 


