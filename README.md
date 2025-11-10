## Authorship
Developed by Erik Kling Åhlén (Student, KTH), 2025
If you use this code in your research, please cite:
> Erik Kling Åhlén (2025). *Thesis-Proteomics-2025: Interactive Analysis Suite for Proteomic Data*. GitHub. https://github.com/ZTKL/Thesis-Proteomics-2025

This repository contains code for the **Interactive Analysis Suite for Proteomic Data**, developed as a part of my Master's thesis *Deciphering the Heterogeneity in Plasma Proteomic Signatures to Reveal Breast Cancer Risk Phenotypes* at KTH in 2025. Earlier development versions were maintained in a private KTH Git repository; this public release contains the finalized scripts prepared for my Master’s thesis.

The suite is designed to provide a interactive workflow facilitating exploratory analysis, visualization and interpretation of proteomic datasets, with a focus on dimensionality reduction and clustering analysis.

## Features:
- **Interactive Workflow for exploratory analysis of Proteomics datasets**
- **Data Quality Control** - tools for handling raw proteomica data and preparing it for analysis.
- **Preprocessing** - Options for filtering, normalisation and adjustments (e.g. residualization).
- **Confounder Identification** - Options for identification of confounding variables including Surrogate Variable Analysis.
- **Stratification** - Identification of dominant trait groups and options for data stratification based on covariate information.
- **Parameter Optimization** - Options for optimization of parameters for included clustering algorithms.
- **Analysis** - Including cluster stability analysis, similarity analysis, correlation analysis, and differential abundance analysis.
- **Visualization** - options for visualization of obtained results.

## Running the Interactive Analysis Suite
Step 1: Load all functions found in the following R scripts to the RStudio environment:
- StudySpecFuncs.R
- HelperFunctions.R
- QCFunctions.R
- HandleNAs.R
- VisualizationFunctions.R
- DimensionReduction.R
- MainRunnerFinal.R
- ConfounderIdentification.R
- StratFuncs.R
- PreProcessing.R
- ClusterParameterOptimization.R
- WrapperFuncs.R
- AnalysisFunctions.R
- DifferentialAbundanceFunctions.R

Step 2: Install and load all packages listed in PackLoader.R 
- Some needs to be installed via BiocManager
- ProtPQN can be found on the SchwenkLab GitHub page
- The correct versions of the packages can be found in the file PackageVersions.txt

Step 3: Prepare data
- npx: protein abundance data (samples as rows, proteins as columns)
  
- sinfo: clinical information (samples as rows)
  -	Required columns:
    o	region (for separating regional cohorts)
    o	x_BC (cancer diagnosis)
    o	sample_id
   	
- binfo: protein information (proteins as rows)
  -	Required columns:
    o	region
    o	panel

Step 4: Start Main Menu script
- run_main_menu(sinfo, binfo, npx)

Notes:
Columns detailing region, x_BC, and panel are used by the setup_variables function found in HelperFunctions.R, this function is applied when entering menus for quality control and preprocessing. If this information is not available in the used dataset it is recommended to create variables select.npx, select.sinfo, and select.binfo duplicating the npx, sinfo and binfo variables. The variables starting with “select.” will be the active variables in the analysis-suite and will remain so by choosing the option “Selected” in source selecting prompts. As the active variables will be changed and overwritten when applying QC, preprocessing and other functions found in the suite, initial creation of the npx, sinfo, and binfo as well as the active variables select.npx, select.sinfo, and select.binfo is recommended. Dimensionality reduction and scaling will create backup variables unreduced.npx and unscaled.npx to be used in downstream analyses, these variables will be overwritten if reapplying such methods.

