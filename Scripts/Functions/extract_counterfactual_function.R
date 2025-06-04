## ---------------------------
##
## Script name: Function to read in counterfactual data
##
## Author: Dr. Joan Dudney
##
## Date Created: 2024-10-22
##
## Copyright (c) Joan Dudney, 2024
## Email: dudney@ucsb.edu
##
## ---------------------------
##
## Notes: Reads in NC files for extraction
##   
##
## ---------------------------

process_nc_file <- function(nc_file, waypoints) {
  # Read the NC file
  r <- terra::rast(nc_file)
  
  # Extract year from filename using regex
  year_range <- stringr::str_extract(nc_file, "\\d{4}_\\d{4}")
  start_year <- as.numeric(stringr::str_extract(year_range, "^\\d{4}"))
  
  # Extract values for waypoints
  extracted_vals <- terra::extract(r, waypoints)
  
  # Create date sequence for the file
  if(start_year == 2011) {
    # Special case for 2011-2019 file
    dates <- seq(as.Date(paste0(start_year, "-01-01")), 
                 as.Date("2019-12-31"), 
                 by = "day")
  } else {
    # Regular case for decade files
    dates <- seq(as.Date(paste0(start_year, "-01-01")), 
                 as.Date(paste0(start_year + 9, "-12-31")), 
                 by = "day")
  }
  
  # Remove the ID column from extracted values
  extracted_vals$ID <- NULL
  
  # Create a data frame with dates
  result_df <- data.frame(
    date = rep(dates, each = nrow(waypoints)),
    point_id = rep(1:nrow(waypoints), length(dates)),
    year = year(dates),
    month = month(dates),
    day = day(dates)
  )
  
  # Add the extracted values
  # Reshape the extracted values to long format
  extracted_long <- tidyr::pivot_longer(extracted_vals, 
                                        cols = everything(),
                                        names_to = "day_number",
                                        values_to = "temperature")
  
  result_df$temperature <- extracted_long$temperature
  
  return(result_df)
}
