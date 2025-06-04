## ---------------------------
##
## Script name: Robustness checks of counterfactual scenarios
##
## Author: Dr. Joan Dudney
##
## Date Created: 2024-08-02
##
## Copyright (c) Joan Dudney, 2024
## Email: dudney@ucsb.edu
##
## ---------------------------
##
## Notes:
##   This script tests whether the main results are robust to different 
##    estimates of the counterfactual temperature scenarios
## ---------------------------


## packages
librarian::shelf(sjPlot, ggeffects, patchwork, tidyverse, broom, progress,
                 lme4, plotrix, ggpubr, mgcv, nlme, fixest, plotrix, egg, ggpmisc,
                 mvtnorm, clubSandwich, rasterVis, broom.mixed, scales, RColorBrewer, splines, zoo)

theme_set(
  theme_bw(base_size = 20)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
# 1. Reading in data and cleaning data
# 
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

select <- dplyr::select

## importing the panel ITRDB data
paneldat_itrdb <- read_csv(here("Data", "SevenSpecies_ITRDB_climatewindows.csv"))

paneldat <- paneldat_itrdb %>%
  mutate(tree_id = paste0(collection_id, "_", tree)) %>% 
  filter(species_id == "pied") %>%
  select(-tree) %>% ## remove duplicated tree ids 
  rename(tree = tree_id, plot = collection_id)

## creating a dataset with tree and plot
treedat <- paneldat %>% 
  select(tree, plot) %>% 
  distinct()

## importing climate data for ITRDB sites
itrdbdat <- read_csv(here("Data", "ITRDB_species_latlon.csv")) %>%
  filter(species_id == "pied") 

## isolating sites for PIED
pied <- itrdbdat %>% 
  select(species_id, collection_id) %>% 
  distinct() %>% 
  rename(plot = collection_id)

## prism data
prism_itrdb <- read_csv(here("Data", "raw_prism_itrdb_sevenDomSp.csv"))

pied_clim <- prism_itrdb %>%
  select(waypoint_id, variable,value, month, year, longitude, latitude) %>% 
  rename(plot = waypoint_id) %>%
  left_join(pied) %>% 
  filter(!is.na(species_id))


climdat <- pied_clim %>% 
  mutate(growing=ifelse(month%in%c(10:12), year+1, year)) %>% 
  filter(!month%in%c(7:9)) %>% ## don't include July, August, September
  filter(growing!=1900) %>% 
  pivot_wider(names_from=variable, values_from=value) %>% 
  group_by(plot, growing) %>% 
  summarize(tmax = mean(tmax, na.rm=T), ppt = sum(ppt, na.rm=T))%>%
  rename(year=growing) %>% 
  ungroup() 


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
# 2. Creating counterfactual data using different detrending approaches
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Base counterfactual function
create_counterfactual <- function(model, data_df, model_name) {
  
  ## Create prediction dataset
  preddat <- data_df
  preddat$predtmax <- predict(model, newdata = preddat)
  
  ## Calculate plot-level mean tmax
  meanpre80 <- preddat %>%
    filter(year<1980) %>% 
    group_by(plot) %>% 
    summarize(meantmax = mean(tmax))
  
  ## Join plot means with main dataset
  mean_counterdat <- preddat %>% 
    left_join(meanpre80)
  
  ## Add year-level heterogeneity to the counterfactual scenario
  ## Formula for each predicted tmax: (tmax - predtmax) + meantmax
  counterdat <- mean_counterdat %>% 
    mutate(tmax_counter = (tmax - predtmax) + meantmax) %>%
    mutate(tmax_counter_historic = ifelse(year > 1979, tmax_counter, tmax),
           model_type = model_name)
  
  return(counterdat)
}

## 1. natural spline (target method)
ns_model <- lm(tmax ~ ns(year, df = 4) + factor(plot), data = climdat)
ns_counter <- create_counterfactual(ns_model, climdat, "Natural Spline")

## 2. GAM model
gam_model <- gam(tmax ~ s(year) + factor(plot), data = climdat)
gam_counter <- create_counterfactual(gam_model, climdat, "GAM")

## 3. linear model
linear_model <- lm(tmax ~ year + factor(plot), data = climdat)
linear_counter <- create_counterfactual(linear_model, climdat, "Linear")

## 4. quadratic model
qu_model <- lm(tmax ~ year + I(year^2) + factor(plot), data = climdat)
qu_counter <- create_counterfactual(qu_model, climdat, "Quadratic")

## 5. smooth spline
spline_fit <- lm(tmax ~ bs(year, knots = c(1980)) + factor(plot), data = climdat)
spline_counter <- create_counterfactual(ns_model, climdat, "Smooth Spline")


## combine all counterfactual datasets
all_counters <- bind_rows(
  ns_counter,
  gam_counter,
  linear_counter,
  qu_counter,
  spline_counter
)

## Visualize the scenarios
yearly_means <- all_counters %>%
  group_by(year, model_type) %>%
  summarize(
    mean_tmax = mean(tmax),
    mean_counter = mean(tmax_counter_historic)
  )

## Plot
ggplot(yearly_means, aes(x = year)) +
  geom_line(aes(y = mean_tmax), color = "black",   size = 1) +
  geom_line(aes(y = mean_counter, color = model_type),  size = 1) +
  geom_vline(xintercept = 1979.5, color = "gray50") +
  geom_smooth(aes(y = mean_counter, color = model_type), se=F)+
  geom_smooth(aes(y = mean_tmax), size = 1)

## Summary statistics (post-1979)
counter_summary <- all_counters %>%
  filter(year > 1979) %>%
  group_by(model_type) %>%
  summarize(
    mean_diff = mean(tmax - tmax_counter_historic),
    max_diff = max(tmax - tmax_counter_historic),
    min_diff = min(tmax - tmax_counter_historic),
    sd_diff = sd(tmax - tmax_counter_historic)
  )

counter_summary

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## creating counterfactual and actual dataframes for prediction
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

actualdat <- climdat %>% 
  left_join(treedat)%>% 
  filter(year<2014)

## GAM dataset
counter_gam <- all_counters %>%
  filter(model_type == "GAM") %>% 
  select(-tmax) %>% 
  rename(tmax = tmax_counter_historic) %>% 
  select(plot, year, ppt, tmax) %>% 
  left_join(treedat)%>% 
  filter(year<2014)

## Linear dataset
counter_linear <- all_counters %>%
  filter(model_type == "Linear") %>% 
  select(-tmax) %>% 
  rename(tmax = tmax_counter_historic) %>% 
  select(plot, year, tmax, ppt)%>% 
  left_join(treedat)%>% 
  filter(year<2014)

## Smooth spline
counter_spline <- all_counters %>%
  filter(model_type == "Smooth Spline") %>% 
  select(-tmax) %>% 
  rename(tmax =  tmax_counter_historic) %>%
  select(plot, year, ppt, tmax)%>%
  left_join(treedat)%>%
  filter(year<2014)

## Quadratic dataset
counter_qu <- all_counters %>%
  filter(model_type == "Quadratic") %>% 
  select(-tmax) %>% 
  rename(tmax = tmax_counter_historic) %>% 
  select(plot, year, ppt, tmax)%>% 
  left_join(treedat)%>% 
  filter(year<2014)

## Target dataset
counter_target <- all_counters %>%
  filter(model_type == "Natural Spline") %>% 
  select(-tmax) %>% 
  rename(tmax = tmax_counter_historic) %>% 
  select(plot, year, ppt, tmax)%>% 
  left_join(treedat)%>% 
  filter(year<2014)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
# 3. MC simulation to compare detrending approaches
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## primary TWFE model
fe_mod <-  feols(rwi ~ tmax * ppt | tree + year,
                 data= paneldat, cluster = ~ plot)

summary(fe_mod)

V_CR1 = vcov(fe_mod)
V_CR1=as.matrix(V_CR1)
coef_vector = fe_mod$coefficients


## randomly drawing coefficients from distribution
draw = rmvnorm(n = 1000, mean = coef_vector, sigma = V_CR1)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MONTE CARLO SIMULATION
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## create an empty dataframe for results
results=data.frame()

pb <- progress_bar$new(
  format = "  Processing [:bar] :percent eta: :eta",
  total = 1000,
  clear = FALSE,
  width = 60
)

for (i in 1:100){
  
  pb$tick()
  
  d <- draw[i,]
  modified_fe_mod <-  fe_mod
  modified_fe_mod$coefficients <- d
  actualdat$vals_actual <-  predict(modified_fe_mod, newdata = actualdat)
  actualdat$vals_gam <-  predict(modified_fe_mod, newdata = counter_gam)
  actualdat$vals_spline <-  predict(modified_fe_mod, newdata = counter_spline)
  actualdat$vals_qu <-  predict(modified_fe_mod, newdata = counter_qu)
  actualdat$vals_linear <-  predict(modified_fe_mod, newdata = counter_linear)
  actualdat$vals_target <-  predict(modified_fe_mod, newdata = counter_target)
  
  results_dat_1980 <- actualdat %>% 
    filter(year >= 1981 & year <= 1991) %>%  
    mutate(tree_diff_gam = vals_actual - vals_gam) %>% 
    mutate(tree_diff_spline = vals_actual - vals_spline) %>% 
    mutate(tree_diff_qu = vals_actual - vals_qu) %>% 
    mutate(tree_diff_linear = vals_actual - vals_linear) %>% 
    mutate(tree_diff_target = vals_actual - vals_target) %>% 
    mutate(iteration = i, period = "1981-1991")
  
  results_dat_2000 <- actualdat %>% 
    filter(year >= 1992 & year <= 2002) %>% 
    mutate(tree_diff_gam = vals_actual - vals_gam) %>% 
    mutate(tree_diff_spline = vals_actual - vals_spline) %>% 
    mutate(tree_diff_qu = vals_actual - vals_qu) %>% 
    mutate(tree_diff_linear = vals_actual - vals_linear) %>%
    mutate(tree_diff_target = vals_actual - vals_target) %>% 
    mutate(iteration = i, period = "1992-2002")
  
  results_dat_2010 <- actualdat %>% 
    filter(year >= 2003 & year <= 2013) %>% 
    mutate(tree_diff_gam = vals_actual - vals_gam) %>% 
    mutate(tree_diff_spline = vals_actual - vals_spline) %>% 
    mutate(tree_diff_qu = vals_actual - vals_qu) %>% 
    mutate(tree_diff_linear = vals_actual - vals_linear) %>%
    mutate(tree_diff_target = vals_actual - vals_target) %>% 
    mutate(iteration = i, period = "2003-2013")
  
  results <- rbind(results, results_dat_1980, results_dat_2000, results_dat_2010)
}

head(results)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# THE FIGURE
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## new dataframe to create figures from
pred_results <- results

## calculate historical mean precip for each plot
tercdat <- climdat %>% 
  filter(year < 1979) %>% 
  group_by(plot) %>% 
  summarize(meanppt = mean(ppt))

tercile_boundaries <- quantile(tercdat$meanppt, probs = c(1/3, 2/3))

tercdat <- tercdat %>%
  mutate(precip_terciles = case_when(
    meanppt < tercile_boundaries[1] ~ "0-33.3%",
    meanppt >= tercile_boundaries[1] & meanppt < tercile_boundaries[2] ~ "33.4-66.7%",
    TRUE ~ "66.8-100%"
  ))

# Sample the dataset to be able to plot it
calc_terc <- pred_results %>%
  left_join(tercdat)


# Calculate the mean and 95% CI for each year and tercile
average_data <- calc_terc %>%
  group_by(period, precip_terciles) %>%
  summarise(
    linear_diff = mean(tree_diff_linear),
    qu_diff = mean(tree_diff_qu),
    spline_diff = mean(tree_diff_spline),
    gam_diff = mean(tree_diff_gam),
    target_diff = mean(tree_diff_target))%>% 
  pivot_longer(
    cols = -c(precip_terciles, period),   
    names_to = "model",                   
    values_to = "mean"                   
  ) %>% 
  mutate(model = str_remove(model, "_diff"))


confdata <- calc_terc %>%
  group_by(period, precip_terciles, iteration) %>%
  summarise(
    mean_diff_linear = mean(tree_diff_linear),
    mean_diff_gam = mean(tree_diff_gam), 
    mean_diff_spline = mean(tree_diff_spline),
    mean_diff_qu = mean(tree_diff_qu),
    mean_diff_target = mean(tree_diff_target) # Changed from mean_diff_targ
  ) %>% 
  group_by(period, precip_terciles) %>% 
  summarise(
    lower_ci_linear = quantile(mean_diff_linear, probs=0.025),
    upper_ci_linear = quantile(mean_diff_linear, probs=0.975),
    lower_ci_gam = quantile(mean_diff_gam, probs=0.025),
    upper_ci_gam = quantile(mean_diff_gam, probs=0.975),
    lower_ci_spline = quantile(mean_diff_spline, probs=0.025),
    upper_ci_spline = quantile(mean_diff_spline, probs=0.975),
    lower_ci_qu = quantile(mean_diff_qu, probs=0.025),
    upper_ci_qu = quantile(mean_diff_qu, probs=0.975),
    lower_ci_target = quantile(mean_diff_target, probs=0.025), # Changed from lower_ci_tar
    upper_ci_target = quantile(mean_diff_target, probs=0.975)  # Changed from upper_ci_tar
  ) %>%
  pivot_longer(
    cols = starts_with("lower_ci_") | starts_with("upper_ci_"),
    names_to = c(".value", "model"),
    names_pattern = "(lower|upper)_ci_(.*)")

results_means_ci <- average_data %>% 
  left_join(confdata) 


## the figure
results_means_ci %>% 
  mutate(model = factor(model, levels = c("target", "linear", "qu", "gam", "spline"))) %>% 
  ggplot(aes(x = factor(period), y = mean, color = model, fill = model)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = lower, ymax = upper), color = "black",
                width = 0.1, 
                position = position_dodge(width = 0.5)) +  # Match dodge width with points
  facet_grid(~precip_terciles)+
  labs(title = "Estimated climate change impacts",
       x = "Time period",
       y = "Δ RWI", color="Method", fill="Method") +
  #scale_color_manual(values = color_palette) +
  #scale_fill_manual(values = color_palette) +
  #guides(fill=F, color=F)+
  scale_color_manual(values = c("#91376F","#303077","#44A894","#EAB94B", "#5e6988"),
                     labels = c("Target", "Linear","Quadratic", "GAM", 
                                "Spline"))+
  scale_fill_manual(values = c("#91376F","#303077","#44A894","#EAB94B", "#5e6988"),
                    labels = c("Target", "Linear","Quadratic", "GAM", 
                               "Spline"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5))



