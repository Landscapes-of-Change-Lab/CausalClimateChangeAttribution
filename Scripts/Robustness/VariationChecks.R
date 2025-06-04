## ---------------------------
##
## Script name: Variation checks
##
## Author: Dr. Joan Dudney
##
## Date Created: 2025-03-29
##
## Copyright (c) Joan Dudney, 2025
## Email: dudney@ucsb.edu
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------


## packages
librarian::shelf(sjPlot, lme4, patchwork, tidyverse, 
                 fixest, clubSandwich, here, lmtest, sandwich, mgcv,
                 marginaleffects)

## setting a figure theme
theme_set(
  theme_bw(base_size = 15)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)

## data to run all analyses
here()
paneldat <- read_csv(here("Data", "cleaned_RWL_climdat.csv")) %>% 
  mutate(ppt = ppt/1000) ## converting mm to m for interpretability


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# EACH MODEL
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mod_fixedeffects = lm(log(value+0.01) ~ factor(year) + factor(tree_id), data= paneldat)
summary(mod_fixedeffects)

mod_year = lm(log(value+0.01) ~ factor(year), data= paneldat)
summary(mod_year)

mod_site = lm(log(value+0.01) ~ factor(tree_id), data= paneldat)
summary(mod_site)

mod_explanatoryvar = lm(log(value+0.01) ~ tmax * ppt, data = paneldat)
summary(mod_explanatoryvar)

mod_all = feols(log(value+0.01) ~ tmax * ppt | tree_id + year, data = paneldat)
summary(mod_all)

precipmod = lm(ppt ~ factor(year), data = paneldat)
summary(precipmod)

precipmod_plot = lm(ppt ~ factor(year) + factor(plot), data = paneldat)
summary(precipmod)

tempmod = lm(tmax ~ factor(year), data = paneldat)
summary(tempmod)

tempmod_plot = lm(tmax ~ factor(year) + factor(plot), data = paneldat)
summary(tempmod_plot)

tempmod_plot = lm(tmax ~  factor(plot), data = paneldat)
summary(tempmod_plot)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# FIGURE of R2 values
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Create a data frame with model names and R² values
r_squared_data <- data.frame(
  model = c("Fixed effects (year + tree_id)", 
            "Year only", 
            "Tree ID only", 
            "Explanatory variables (tmax * ppt)",
            "Temperature ~ year",
            "Temperature ~ year + plot",
            "Precipitation ~ year",
            "Precipitation ~ year + plot"),
  adj_r_squared = c(
    summary(mod_fixedeffects)$adj.r.squared,
    summary(mod_year)$adj.r.squared,
    summary(mod_site)$adj.r.squared,
    summary(mod_explanatoryvar)$adj.r.squared,
    summary(tempmod)$adj.r.squared,
    summary(tempmod_plot)$adj.r.squared,
    summary(precipmod)$adj.r.squared,
    summary(precipmod_plot)$adj.r.squared
  )
)

# Sort data by adjusted R² values 
r_squared_data <- r_squared_data[order(-r_squared_data$adj_r_squared),]


# Create the dot plot adjusted R² 
ggplot(r_squared_data, aes(x = adj_r_squared, y = reorder(model, adj_r_squared))) +
  geom_segment(aes(x = 0, xend = adj_r_squared, yend = model), 
               color = "#a6bddb", size = 0.75) +  
  geom_point(size = 4, color = "#2c7fb8") +  
  labs(
       x = "Adjusted R-squared Value",
       y = "") +  
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 19),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
  ) +
  geom_text(aes(label = sprintf("%.3f", adj_r_squared)), 
            hjust = -0.3, size = 4.5) +  
  xlim(0, max(r_squared_data$adj_r_squared) * 1.15)  




