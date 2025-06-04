## ---------------------------
##
## Script name: Building the pec chart 
##
## Author: Dr. Joan Dudney
##
## Date Created: 2024-09-15
##
## Copyright (c) Joan Dudney, 2024
## Email: dudney@ucsb.edu
##
## ---------------------------
##
## Notes: This code estimates different models with 
##    specifications to build the specification chart 
##    which is Figure 7 in the manuscript

## This code was originally developed by Ariel Ortiz Bobea: 
## https://github.com/ArielOrtizBobea/spec_chart
## ---------------------------



## packages
librarian::shelf(sjPlot, lme4, patchwork, tidyverse, 
                 fixest, clubSandwich, here, lmtest, sandwich, mgcv,
                 marginaleffects)

select <- dplyr::select

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# A. Reading in data and creating the panel dataset

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## source function required for this analysis
source(here("Scripts", "Functions", "spec_chart_function.R"))

## itrdb data
paneldat_itrdb <- read_csv(here("Data", "SevenSpecies_ITRDB_climatewindows.csv")) %>% 
  filter(ppt<1000) %>% 
  mutate(ppt = ppt/1000, pptSummer = pptSummer/1000, 
         ppt_an = ppt_an/1000, laggedprecip = laggedprecip/1000) ## converting mm to m; more interpretable coefs

paneldat <- paneldat_itrdb %>%
  mutate(tree_id = paste0(collection_id, "_", tree)) %>% 
  filter(species_id == "pied") %>%
  select(-tree) %>% ## remove duplicated tree ids 
  rename(tree = tree_id, plot = collection_id)

length(unique(paneldat$tree))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# B. Create function to extract coefficients and empty dataframe to bind model results

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Function to extract coefficients and SE
extract_coef_se <- function(model) {
  slopes_df <- avg_slopes(model) %>% 
    as_tibble() %>% 
    print()
  slopes_df %>% 
    filter(term %in% c("tmax", "tmaxSummer", "tmax_an")) %>% 
    select(term, estimate, std.error)
}

## empty dataframe
specs <- tibble(coef=NaN, 
                    se=NaN,
                    clim_lin = NaN,
                    clim_quad = NaN,
                    clim_int = NaN,
                    clim_lag = FALSE,
                    ww_wy = NaN,
                    ww_an = NaN,
                    ww_sum = NaN,
                    struc_yearfe = NaN,
                    struc_plotfe = NaN,
                    struc_treefe = NaN,
                    struc_ri = NaN,
                    struc_rs = NaN,
                    drop_max = NaN,
                    drop_min = NaN,
                    stder_clust = NaN,
                    stder_het = NaN,
                    stder_conley = NaN) 


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# C. Run each different model and extract coefficients for spec chart

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1. The target panel model 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fe_mod <-  feols(rwi ~ tmax * ppt | tree + year,
                 data= paneldat, cluster = ~ plot)

summary(fe_mod)


## extract coefficients
fe_coef <- extract_coef_se(fe_mod)

## add correct model specification
new_row <-  tibble(coef=fe_coef$estimate, 
                       se=fe_coef$std.error,
                       clim_lin = TRUE,
                       clim_quad = FALSE,
                       clim_int = TRUE,
                       clim_lag = FALSE,
                       ww_wy = TRUE,
                       ww_an = FALSE,
                       ww_sum = FALSE,
                       struc_yearfe = TRUE,
                       struc_treefe = TRUE,
                       struc_plotfe = FALSE,
                       struc_ri = FALSE,
                       struc_rs = FALSE,
                       drop_max = FALSE,
                       drop_min = FALSE,
                       stder_clust = TRUE,
                       stder_het = FALSE,
                       stder_conley = FALSE)

## combine data
specs <- rbind(specs, new_row)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2. Conley standard errors (spatial autocorrelation)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe_mod_con <-  feols(rwi ~ tmax * ppt | tree + year,
                 data= paneldat, vcov = "conley")

## extract coefficients
fe_coef_con <- extract_coef_se(fe_mod_con)

## add correct model specification
new_row <-  tibble(coef=fe_coef_con$estimate, 
                       se=fe_coef_con$std.error,
                       clim_lin = TRUE,
                       clim_quad = FALSE,
                       clim_int = TRUE,
                       clim_lag = FALSE,
                       ww_wy = TRUE,
                       ww_an = FALSE,
                       ww_sum = FALSE,
                       struc_yearfe = TRUE,
                       struc_treefe = TRUE,
                       struc_plotfe = FALSE,
                       struc_ri = FALSE,
                       struc_rs = FALSE,
                       drop_max = FALSE,
                       drop_min = FALSE,
                       stder_clust = FALSE,
                       stder_het = FALSE,
                       stder_conley = TRUE)

## combine data
specs <- rbind(specs, new_row)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3. heteroskedasticity-robust standard errors
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe_mod_hetero <-  feols(rwi ~ tmax * ppt | tree + year,
                        data= paneldat, vcov = "hetero")

## extract coefficients
fe_hetero_coef <- extract_coef_se(fe_mod_hetero)

## add correct model specification
new_row <-  tibble(coef=fe_hetero_coef$estimate, 
                       se=fe_hetero_coef$std.error,
                       clim_lin = TRUE,
                       clim_quad = FALSE,
                       clim_int = TRUE,
                       clim_lag = FALSE,
                       ww_wy = TRUE,
                       ww_an = FALSE,
                       ww_sum = FALSE,
                       struc_yearfe = TRUE,
                       struc_treefe = TRUE,
                       struc_plotfe = FALSE,
                       struc_ri = FALSE,
                       struc_rs = FALSE,
                       drop_max = FALSE,
                       drop_min = FALSE,
                       stder_clust = FALSE,
                       stder_het = TRUE,
                       stder_conley = FALSE)

## combine data
specs <- rbind(specs, new_row)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3. Panel model with quadratic terms
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe_mod_quad <-  feols(rwi ~ tmax + ppt + I(tmax^2) + I(ppt^2) | tree + year,
                      data= paneldat, cluster = ~ plot)
#summary(fe_mod_quad)

## extract coefficients
fe_quad_coef <- extract_coef_se(fe_mod_quad)

## add correct model specification
new_row <-  tibble(coef=fe_quad_coef$estimate, 
                       se=fe_quad_coef$std.error,
                       clim_lin = FALSE,
                       clim_quad = TRUE,
                       clim_int = FALSE,
                       clim_lag = FALSE,
                       ww_wy = TRUE,
                       ww_an = FALSE,
                       ww_sum = FALSE,
                       struc_yearfe = TRUE,
                       struc_treefe = TRUE,
                       struc_plotfe = FALSE,
                       struc_ri = FALSE,
                       struc_rs = FALSE,
                       drop_max = FALSE,
                       drop_min = FALSE,
                       stder_clust = TRUE,
                       stder_het = FALSE,
                       stder_conley = FALSE)

## combine data
specs <- rbind(specs, new_row)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4. Panel model without interaction effects 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fe_mod_no_int<-  feols(rwi ~ tmax + ppt | tree + year,
                       data= paneldat, cluster = ~ plot)

#summary(fe_mod_no_int)

## extract coefficients
fe_no_int_coef <- extract_coef_se(fe_mod_no_int)


## add correct model specification
new_row <-  tibble(coef=fe_no_int_coef$estimate, 
                       se=fe_no_int_coef$std.error,
                       clim_lin = TRUE,
                       clim_quad = FALSE,
                       clim_int = FALSE,
                       clim_lag = FALSE,
                       ww_wy = TRUE,
                       ww_an = FALSE,
                       ww_sum = FALSE,
                       struc_yearfe = TRUE,
                       struc_treefe = TRUE,
                       struc_plotfe = FALSE,
                       struc_ri = FALSE,
                       struc_rs = FALSE,
                       drop_max = FALSE,
                       drop_min = FALSE,
                       stder_clust = TRUE,
                       stder_het = FALSE,
                       stder_conley = FALSE)

## add data
specs <- rbind(specs, new_row)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5. Panel model with lagged effects
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fe_mod_lag <-  feols(rwi ~ tmax * ppt + laggedtmax + laggedprecip | tree + year,
                     data= paneldat, cluster = ~ plot)

summary(fe_mod_lag)

## extact coefficients
fe_lag_coef <- extract_coef_se(fe_mod_lag)

## add correct model specification
new_row <-  tibble(coef=fe_lag_coef$estimate, 
                       se=fe_lag_coef$std.error,
                       clim_lin = TRUE,
                       clim_quad = FALSE,
                       clim_int = TRUE,
                       clim_lag = TRUE,
                       ww_wy = TRUE,
                       ww_an = FALSE,
                       ww_sum = FALSE,
                       struc_yearfe = TRUE,
                       struc_treefe = TRUE,
                       struc_plotfe = FALSE,
                       struc_ri = FALSE,
                       struc_rs = FALSE,
                       drop_max = FALSE,
                       drop_min = FALSE,
                       stder_clust = TRUE,
                       stder_het = FALSE,
                       stder_conley = FALSE)

## combine data
specs <- rbind(specs, new_row)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 6. Panel model with annual weather
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe_mod_an <-  feols(rwi ~ tmax_an * ppt_an | tree + year,
                    data= paneldat, cluster = ~ plot)


## extract coefficients
fe_an_coef <- extract_coef_se(fe_mod_an)

## add correct model specification
new_row <-  tibble(coef=fe_an_coef$estimate, 
                       se=fe_an_coef$std.error,
                       clim_lin = TRUE,
                       clim_quad = FALSE,
                       clim_int = TRUE,
                       clim_lag = FALSE,
                       ww_wy = FALSE,
                       ww_an = TRUE,
                       ww_sum = FALSE,
                       struc_yearfe = TRUE,
                       struc_treefe = TRUE,
                       struc_plotfe = FALSE,
                       struc_ri = FALSE,
                       struc_rs = FALSE,
                       drop_max = FALSE,
                       drop_min = FALSE,
                       stder_clust = TRUE,
                       stder_het = FALSE,
                       stder_conley = FALSE)

## combine data
specs <- rbind(specs, new_row)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 7. Panel model with summer weather
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe_mod_summer <-  feols(rwi ~ tmaxSummer * pptSummer | tree + year,
                        data= paneldat, cluster = ~ plot)


## extract coefficients
fe_summer_coef <- extract_coef_se(fe_mod_summer)

## add correct model specification
new_row <-  tibble(coef=fe_summer_coef$estimate, 
                       se=fe_summer_coef$std.error,
                       clim_lin = TRUE,
                       clim_quad = FALSE,
                       clim_int = TRUE,
                       clim_lag = FALSE,
                       ww_wy = FALSE,
                       ww_an = FALSE,
                       ww_sum = TRUE,
                       struc_yearfe = TRUE,
                       struc_treefe = TRUE,
                       struc_plotfe = FALSE,
                       struc_ri = FALSE,
                       struc_rs = FALSE,
                       drop_max = FALSE,
                       drop_min = FALSE,
                       stder_clust = TRUE,
                       stder_het = FALSE,
                       stder_conley = FALSE)

## add data
specs <- rbind(specs, new_row)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 8. LM model
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotfe_mod <- feols(rwi ~ tmax * ppt | tree + year,
                    data = paneldat, cluster = ~ plot)

## extract coefficients
plotfe_coef <- extract_coef_se(plotfe_mod)

## add correct model specification
new_row <-  tibble(coef=plotfe_coef$estimate, 
                       se=plotfe_coef$std.error,
                       clim_lin = TRUE,
                       clim_quad = FALSE,
                       clim_int = TRUE,
                       clim_lag = FALSE,
                       ww_wy = TRUE,
                       ww_an = FALSE,
                       ww_sum = FALSE,
                       struc_yearfe = TRUE,
                       struc_treefe = FALSE,
                       struc_plotfe = TRUE,
                       struc_ri = FALSE,
                       struc_rs = FALSE,
                       drop_max = FALSE,
                       drop_min = FALSE,
                       stder_clust = TRUE,
                       stder_het = FALSE,
                       stder_conley = FALSE)

## add data
specs <- rbind(specs, new_row)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 9. LMM with random intercepts
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rei_mod = lmer(rwi ~ tmax *  ppt +  (1|tree/plot), 
               control = lmerControl(optimizer = "bobyqa",
                                     optCtrl = list(maxfun = 100000)), 
               data = paneldat)

#tab_model(rei_mod)

## extract coefficients
rei_coef <- extract_coef_se(rei_mod)


## add correct model specification
new_row <-  tibble(coef=rei_coef$estimate, 
                       se=rei_coef$std.error,
                       clim_lin = TRUE,
                       clim_quad = FALSE,
                       clim_int = TRUE,
                       clim_lag = FALSE,
                       ww_wy = TRUE,
                       ww_an = FALSE,
                       ww_sum = FALSE,
                       struc_yearfe = FALSE,
                       struc_treefe = FALSE,
                       struc_plotfe = FALSE,
                       struc_ri = TRUE,
                       struc_rs = FALSE,
                       drop_max = FALSE,
                       drop_min = FALSE,
                       stder_clust = FALSE, 
                       stder_het = FALSE,
                       stder_conley = FALSE)

## combine data
specs <- rbind(specs, new_row)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 10. Just random slopes
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Updated model with scaled variables and increased iterations
res_mod <- lmer(rwi ~ tmax * ppt + 
                  (0 + tmax | tree) +
                  (0 + ppt | tree), 
                control = lmerControl(optimizer = "bobyqa",
                                      optCtrl = list(maxfun = 100000)), 
                data = paneldat)

## extract coefficients
res_coef <- extract_coef_se(res_mod)

## add correct model specification
new_row <-  tibble(coef=res_coef$estimate, 
                       se=res_coef$std.error,
                       clim_lin = TRUE,
                       clim_quad = FALSE,
                       clim_int = TRUE,
                       clim_lag = FALSE,
                       ww_wy = TRUE,
                       ww_an = FALSE,
                       ww_sum = FALSE,
                       struc_yearfe = FALSE,
                       struc_treefe = FALSE,
                       struc_plotfe = FALSE,
                       struc_ri = FALSE, 
                       struc_rs = TRUE,
                       drop_max = FALSE,
                       drop_min = FALSE,
                       stder_clust = FALSE,
                       stder_het = FALSE,
                       stder_conley = FALSE)

## combine data
specs <- rbind(specs, new_row)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 11. Panel model without year effects
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe_mod_noyr <-  feols(rwi ~ tmax * ppt | tree,
                      data= paneldat, cluster = ~ plot)

#summary(fe_mod_noyr)

## extract coefficients
fe_noyr_coef <- extract_coef_se(fe_mod_noyr)

## add correct model specification
new_row <-  tibble(coef=fe_noyr_coef$estimate, 
                       se=fe_noyr_coef$std.error,
                       clim_lin = TRUE,
                       clim_quad = FALSE,
                       clim_int = TRUE,
                       clim_lag = FALSE,
                       ww_wy = TRUE,
                       ww_an = FALSE,
                       ww_sum = FALSE,
                       struc_yearfe = FALSE,
                       struc_treefe = TRUE,
                       struc_plotfe = FALSE,
                       struc_ri = FALSE,
                       struc_rs = FALSE,
                       drop_max = FALSE,
                       drop_min = FALSE,
                       stder_clust = TRUE,
                       stder_het = FALSE,
                       stder_conley = FALSE)

## combined ata
specs <- rbind(specs, new_row)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Jackknife --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## pull plot values
plots <- paneldat %>% pull(plot) %>% unique()

## empty dataframe for for-loop values
drop_coefs <- tibble()

## for loop
for (p in plots) {
  drop_dat <- paneldat %>% filter(plot != p)
  
  drop_mod <-  feols(rwi ~ tmax * ppt | tree + year,
                     data= drop_dat, cluster = ~ plot)
  drop_coef <- extract_coef_se(drop_mod)
  drop_coefs <- rbind(drop_coefs, drop_coef)
}

## extract min and max values 
min_mod <- drop_coefs %>%  slice(which.min(estimate))
max_mod <- drop_coefs %>%  slice(which.max(estimate))

## add correct model specification
## max value
new_row <-  tibble(coef=max_mod$estimate, 
                       se=max_mod$std.error,
                       clim_lin = TRUE,
                       clim_quad = FALSE,
                       clim_int = TRUE,
                       clim_lag = FALSE,
                       ww_wy = TRUE,
                       ww_an = FALSE,
                       ww_sum = FALSE,
                       struc_yearfe = TRUE,
                       struc_treefe = TRUE,
                       struc_plotfe = FALSE,
                       struc_ri = FALSE,
                       struc_rs = FALSE,
                       drop_max = TRUE,
                       drop_min = FALSE,
                       stder_clust = TRUE,
                       stder_het = FALSE,
                       stder_conley = FALSE)

## combine data
specs <- rbind(specs, new_row)

## min value
new_row <-  tibble(coef=min_mod$estimate, 
                       se=min_mod$std.error,
                       clim_lin = TRUE,
                       clim_quad = FALSE,
                       clim_int = TRUE,
                       clim_lag = FALSE,
                       ww_wy = TRUE,
                       ww_an = FALSE,
                       ww_sum = FALSE,
                       struc_yearfe = TRUE,
                       struc_treefe = TRUE,
                       struc_plotfe = FALSE,
                       struc_ri = FALSE,
                       struc_rs = FALSE,
                       drop_max = FALSE,
                       drop_min = TRUE,
                       stder_clust = TRUE,
                       stder_het = FALSE,
                       stder_conley = FALSE)

specs <- rbind(specs, new_row)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create spec chart --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Create figure
highlight_n <- 1

## Define label structure
labels <- list("Climate relationships" = c("Linear", "Quadratic", "Interaction", "Lagged"),
               "Weather window" = c("Water year", "Annual", "Summer"),
               "Model structure" = c("Year FEs", "Tree FEs", "Plot FEs", "Random intercepts", "Random slopes"),
               "Jackknife by plot" = c("Minimum effect", "Maximum effect"),
               "Standard errors" = c("Clustered at plot", "Heteroskedasticity robust", "Conley"))

specs <- specs %>% drop_na()

specs_df <- as.data.frame(specs)

nrow(specs_df)
sum(c(4, 3, 5, 2, 3))

par(oma=c(1,0,1,1))

robustness_fig <- schart(specs_df, labels, highlight=highlight_n, order = "asis", 
                         # cex=(1.2),fonts=c(2,3),
                         # heights = c(.4,.6),
                         # n=c(13), ci = c(.95),
                         n = c(3, 3, 2, 4, 2), ci = c(.95),
                         ylab = "Average marginal effect of\ntemperature on growth",
                         col.est = c("grey80", "dodgerblue4"),
                         col.dot = c("grey60","grey95","grey95","dodgerblue4"),
                         bg.dot = c("grey60","grey95","grey95","dodgerblue4"),
                         lwd.symbol = 1)


text(x = 1, y = -0.007, label = "A",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 2, y = -0.007, label = "B",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 3, y = -0.007, label = "C",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 5, y = -0.007, label = "D",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 6, y = -0.007, label = "E",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 7, y = -0.007, label = "F",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 9, y = -0.007, label = "G",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 10, y = -0.007, label = "H",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 12, y = -0.007, label = "I",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 13, y = -0.007, label = "J",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 14, y = -0.007, label = "K",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 15, y = -0.007, label = "L",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 17, y = -0.007, label = "M",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 18, y = -0.007, label = "N",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

#dev.off()




