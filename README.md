# Causal Climate Change Attribution
Replication data and code for the case study described in Dudney et al., 2025, "A causal inference framework for climate change attribution in ecology"
DOI: https://doi.org/10.5281/zenodo.15611152

## Instructions for replication
To replicate the analysis, readers should complete the following steps:
- Clone repository: This repository contains all code and data needed for replication.
- Create an R project: Open RStudio and create a new project in the directory where you cloned the repository.
- Set up computing environment: Run renv::restore() in the R console to install all required packages and set up the computing environment.
- Run the scripts: The "/Scripts/full_analysis.R" script runs each of the analysis scripts in the correct order. 
Note that the full analysis can take several hours to complete on typical personal computers. Decrease the number of iterations (n_mc) to run more quickly.
Individual scripts are stored in the "Scripts/"directory.
- Explore the outputs: All output figures and tables will be saved to the "Output/" directory. 
