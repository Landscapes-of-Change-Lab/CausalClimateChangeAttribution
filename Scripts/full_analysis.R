#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##
## Script name: Full analysis
##
## Author: Dr. Joan Dudney and Robert Heilmayr
##
## Date Created: 2025-6-5
##
## Copyright (c) Joan Dudney, 2024
## Email: dudney@ucsb.edu
##
## Notes:
##   This script runs the full workflow associated with the case study in 
##    Dudney et al, 2025
##
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ------ Set monte carlo simulation parameters  ------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Set number of monte carlo simulations to run
n_mc = 1000  # Paper uses 1000, but this can take several hours to run.

# Set random seed for reproducibility
random_seed = 93105
set.seed(random_seed)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ------ Create directory to store outputs  ------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dir.create(here("Output"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ------ Run primary analysis  ------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run primary analysis (Figure 6) and test robustness against precip changes (Figure S2)
source("Scripts/Primary analysis/CounterfactualclimatechangeITRDB.R")

# Generate marginal effects plot (Figure 5, Table S2)
source("Scripts/Primary analysis/MarginalEffectsPlot.R")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ------ Run robustness tests  ------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Test robustness against alternate detrending approaches to counterfactual climate (Figure S1)
source("Scripts/Robustness/Robustness_detrending.R")

# Test robustness across multiple model specifications (Figure 7)
source("Scripts/Robustness/Robustness_specChart.R")

# Test for spatial autocorrelation (Table S3) (Warning: can take >1 hour to run)
source("Scripts/Robustness/SpatialAutocorrelationTest.R")
