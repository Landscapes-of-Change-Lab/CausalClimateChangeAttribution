## ---------------------------
##
## Script name: Climate trend tests
##
## Author: Dr. Joan Dudney
##
## Date Created: 2024-09-07
##
## Copyright (c) Joan Dudney, 2024
## Email: dudney@ucsb.edu
##
## ---------------------------
##
## Notes: Tests for different linear and spline trends
##   in the climate data
##
## ---------------------------
## packages
librarian::shelf(sjPlot, ggeffects, patchwork, tidyverse, 
                 lme4, plotrix, ggpubr, mgcv, nlme, fixest, plotrix, egg, ggpmisc,
                 mvtnorm, clubSandwich, rasterVis, broom.mixed, scales,RColorBrewer, splines, zoo, here)

select <- dplyr::select


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Importing and cleaning data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## importing ITRDB data
itrdbdat <- read_csv(here("Data", "ITRDB_species_latlon.csv")) %>%
  filter(species_id == "pied") 

## isolating sites for PIED
pied <- itrdbdat %>% 
  select(species_id, collection_id) %>% 
  distinct()

## prism data
prism_itrdb <- read_csv(here("Data", "raw_prism_itrdb_sevenDomSp.csv"))

pied_clim <- prism_itrdb %>%
  select(waypoint_id, variable,value, month, year, longitude, latitude) %>% 
  rename(collection_id = waypoint_id) %>%
  left_join(pied) %>% 
  filter(!is.na(species_id))

## all climate data
climdat <- pied_clim %>% 
  mutate(growing=ifelse(month%in%c(10:12), year+1, year)) %>% 
  filter(!month%in%c(7:9)) %>% ## don't include July, August, September
  filter(growing!=1900) %>% 
  pivot_wider(names_from=variable, values_from=value) %>% 
  group_by(collection_id, growing) %>% 
  summarize(tmax = mean(tmax, na.rm=T), ppt = sum(ppt, na.rm=T))%>%
  rename(year=growing) %>% 
  ungroup() %>% 
  rename(plot = collection_id)

## visualize
climdat %>% 
  ggplot(aes(x=year, y=tmax))+
  geom_smooth()

## climate data
climdat <- climdat %>% 
  filter(year<2014) %>% 
  filter(year>1980)

## model
feols(tmax ~ year | plot, data = climdat)

feols(ppt ~ year | plot, data = climdat)
