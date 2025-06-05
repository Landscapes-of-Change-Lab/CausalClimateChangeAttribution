# CausalClimateChangeAttribution
Replication data and code for the case study described in Dudney et al., 2025, "A causal inference framework for climate change attribution in ecology"

# Instructions for replication
To replicate the analysis, readers should complete the following steps:
- Clone repository: This repository contains all code and data needed for replication.
- Create an R project: Open RStudio and create a new project in the directory where you cloned the repository.
### Question for Joan: Are you open to using renv for package management
- Install librarian: The scripts in this repository use librarian to manage package dependencies. Run install("librarian") in an R console within your project.
- Run the scripts: The "/Scripts/full_analysis.R" script runs each of the analysis scripts in the correct order. 
Note that this full analysis can take several hours to complete on typical personal computers. 
Individual scripts are stored in the "/Scripts/"directory.
- Explore the outputs: All output figures will be saved to the "/Output/Figures/" directory. 
### Question for Joan: Do you agree that we should be saving scripts out within the directory?