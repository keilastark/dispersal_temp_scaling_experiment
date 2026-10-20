This directory contains the code and data needed to reproduce the results from the manuscript: Dispersal modifies temperature dependence of ecosystem
productivity after an extreme thermal fluctuation. 
We ran an outdoor aquatic mesocosm experiment in the summer of 2018 to test the hypothesis that different dispersal magnitudes modify the thermal sensitivity of Gross Primary Productivity by redistributing temperature-induced differences in biomass, community size spectra, and thermal phenotypes. 


**"Data" folder**

This folder contains data used in the manuscript. Some data files are explicitly called in the R scripts, while others contain extra supplementary data. The most important files are:

"metacom_exp_data_final.csv": Master datasheet for autotroph and zooplankton biomass, GPP, total Nitrogen, and phytoplankton community size structure. For GLMMs/SEM. 
"phytoplankton_diversity_data.csv": Phytoplankton community data (estimated with visual key-based ID of flow cytometry images)
"Bacterial_16S_data.csv": Bacterial community data (estimated with 16S)
"Zoop_data.csv": Zooplankton community data (estimated by microscopy)


**"Scripts" folder**

R scripts are roughly organized by analysis and data type and in the order presented in the manuscript. 



Hypothesis1_LMMS.R runs week-specific Bayesian GLMMs to test the hypothesis that the slope describing the temperature vs GPP relationship is modified by dispersal. The model is described in the Methods section of the manuscript.

Hypothesis2_SEM_and_DAG.R runs week-specific Structural Equation Modelling, and reproduces the Directed Acyclic Graphs in Fig. S9.

Hypothesis2_analyses.R runs the analyses and produces figures for testing Hypothesis 2a-c. These include LMMs, Levene's tests, multivariate PERMDISP tests of treatment effects on community biomass, size structure, taxonomic composition, and thermal phenotype composition.

Supplementary_Analyses.R includes analyses and figures in the Supporting Information, including the reaeration correction, sensitivity analyses for LMMs used as main hypothesis tests, and supplementary data visualizations. 

processing_flux_data.R is a supplementary file showing how GPP was calculated from raw oxygen measurements. 
