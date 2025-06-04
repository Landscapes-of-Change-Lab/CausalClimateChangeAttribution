## ---------------------------
##
## Script name: Moran's I Test
##
## Author: Dr. Joan Dudney
##
## Date Created: 2024-09-11
##
## Copyright (c) Joan Dudney, 2024
## Email: dudney@ucsb.edu
##
## ---------------------------
##
## Notes: This code creates a function that tests for 
##  spatial autocorrelation of residuals in our fixest 
##  panel model of RWI
##
## ---------------------------

librarian::shelf(fixest, spdep, sf,units, here)


paneldat_itrdb <- read_csv(here("Data", "SevenSpecies_ITRDB_climatewindows.csv")) %>% 
  filter(ppt<1000) %>% 
  mutate(ppt = ppt/1000, pptSummer = pptSummer/1000, 
         ppt_an = ppt_an/1000, laggedprecip = laggedprecip/1000) ## converting mm to m; more interpretable coefs

paneldat <- paneldat_itrdb %>%
  mutate(tree_id = paste0(collection_id, "_", tree)) %>% 
  filter(species_id == "pied") %>%
  select(-tree) %>% ## remove duplicated tree ids 
  rename(tree = tree_id, plot = collection_id) %>% 
  rename(lat = latitude, long = longitude)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# The panel model 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe_mod <-  feols(rwi ~ tmax * ppt | tree + year,
                 data= paneldat, cluster = ~ plot)

summary(fe_mod)


# Function that test for spatial autocorrelation in fixest residuals
test_spatial_autocorrelation <- function(model, data, lat = "lat", lon = "long", dist_km = 20) {
  
  # Extract residuals
  residuals <- residuals(model)
  
  # Convert data to sf object
  coords_sf <- st_as_sf(data, coords = c(lon, lat), crs = 4326)  # WGS84
  
  # Transform to a projected CRS for distance calculation
  coords_proj <- st_transform(coords_sf, 3857)  # Web Mercator projection
  
  # Create a spatial weights matrix
  
  cat("Creating spatial weights matrix...\n")
  coords_matrix <- st_coordinates(coords_proj)
  
  # Use dnearneigh to create neighbor list
  nb <- spdep::dnearneigh(coords_matrix, 0, dist_km * 1000)  # Convert km to meters
  
  # Diagnostic prints
  cat("Structure of nb:\n")
  print(str(nb))
  cat("Length of nb:", length(nb), "\n")
  cat("Class of nb:", class(nb), "\n")
  
  # Convert to listw object with error handling
  cat("Converting to listw object...\n")
  tryCatch({
    lw <- spdep::nb2listw(nb, style = "W")
  }, error = function(e) {
    cat("Error in nb2listw:", conditionMessage(e), "\n")
    cat("Attempting to fix nb...\n")
    nb <- spdep::make.sym.nb(nb)
    lw <- spdep::nb2listw(nb, style = "W")
  })
  
  # Perform Moran's I test
  cat("Performing Moran's I test...\n")
  moran_test <- moran.test(residuals, lw)
  
  return(moran_test)
}


# Fit the fixest model

cat("Fitting the model...\n")

#model <- feols(value ~ tmax + ppt + year | plot_id_needle, data = rwidat)
model = fe_mod

# Test for spatial autocorrelation within 20 km

cat("Testing for spatial autocorrelation...\n")
#autocorrelation_result <- test_spatial_autocorrelation(model, rwidat, lat = "lat", lon = "long", dist_km = 20)

autocorrelation_result <- test_spatial_autocorrelation(model, paneldat, lat = "lat", lon = "long", dist_km = 20)

# Print the result
print(autocorrelation_result)


