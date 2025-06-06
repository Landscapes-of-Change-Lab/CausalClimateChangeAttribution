#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##
## Script name: Marginal effects of climate variables
##
## Author: Dr. Joan Dudney
##
## Date Created: 2024-10-20
##
## Copyright (c) Joan Dudney, 2024
## Email: dudney@ucsb.edu
##
##
## Notes: Estimates a fixed effects panel model of PIED growth and
##   plots the marginal effects of temp and precip 
##
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Packages and set-up ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

librarian::shelf(sjPlot, lme4, patchwork, tidyverse, ggeffects,
                 fixest, clubSandwich, here, lmtest, sandwich, marginaleffects)

theme_set(
  theme_bw(base_size = 15)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)

select <- dplyr::select


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Data cleaning ------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## ITRDB data
paneldat <- read_csv(here("Data", "paneldata_ITRDB_PIED.csv")) %>% 
  mutate(ppt = ppt/1000, pptSummer = pptSummer/1000, 
         ppt_an = ppt_an/1000, laggedprecip = laggedprecip/1000) ## converting mm to m; more interpretable coefs


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# The panel model ------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe_mod <-  feols(rwi ~ tmax * ppt | tree + year,
                 data= paneldat, cluster = ~ plot)

summary(fe_mod)


tab_model(fe_mod, digits = 4, show.ci = F, show.se = T)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Margins plot -------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Compute temperature marginal effects
mtemp <- slopes(fe_mod, 
              variables = "tmax",
              by = "ppt",
              newdata = datagrid(ppt = seq(min(paneldat$ppt), max(paneldat$ppt), length.out = 100)))

# create the marginal effects plot of temperature
marg_temp = ggplot(mtemp, aes(x = ppt*1000, y = estimate, color)) +
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_line(color="#91376F") +
  #ylim(-.03, .08)+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color="#91376F", fill="#91376F") +
  labs(x = "Precipitation (mm)", 
       y = "Marginal effect of temperature")+
  scale_x_continuous(breaks = seq(0, 3800, by = 500))


# Compute marginal effects of precip
mprecip <- slopes(fe_mod, 
                variables = "ppt",
                by = "tmax",
                newdata = datagrid(tmax = seq(min(paneldat$tmax), max(paneldat$tmax), length.out = 100)))

# The precipitation plot
marg_precip = ggplot(mprecip, aes(x = tmax, y = estimate)) +
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_line(color= "#303077") +
  #ylim(-.2, .45)+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color="#303077", fill = "#303077") +
  labs(x = "Temperature (°C)", 
       y = "Marginal effect of precipitation")


## visualizing the two plots together
marg_temp + marg_precip +   plot_annotation(tag_levels = "A")
 



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Adding density distribution panels to both plots ------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Create a precipitation density plot
ppt_density <- ggplot(paneldat, aes(x = ppt*1000)) +
  geom_density(fill = "#91376F", alpha = 0.6) +
  scale_x_continuous(breaks = seq(0, 3800, by = 500)) +
  scale_y_continuous(
    breaks = function(x) pretty(x, n = 3),
    labels = function(x) format(x, scientific = FALSE)  # Add this line
  ) +
  ggtitle("A") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    plot.margin = margin(5, 0, -5, 0),  # Negative bottom margin to reduce gap
    panel.border = element_blank(),
    axis.line.x = element_line(color = "gray70"),
    axis.line.y = element_line(color = "gray70"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0, vjust = 1),
    plot.title.position = "plot"
  ) +
  labs(y = "Density")


# Combine with temperature marginal effects plot
combined_temp_plot <- ppt_density / marg_temp +
  plot_layout(heights = c(2, 5), guides = "collect") +
  plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5))) &
  theme(plot.margin = margin(5, 5, 5, 5))


# Create a density plot for temperature
tmax_density <- ggplot(paneldat, aes(x = tmax)) +
  geom_density(fill = "#303077", alpha = 0.6) +
  scale_y_continuous(breaks = function(x) pretty(x, n = 3)) +
  ggtitle("B") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    plot.margin = margin(5, 0, -5, 0),  # Negative bottom margin to reduce gap
    panel.border = element_blank(),
    axis.line.x = element_line(color = "gray70"),
    axis.line.y = element_line(color = "gray70"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0, vjust = 1),
    plot.title.position = "plot"
  ) +
  labs(y = "Density")


# Combine with precipitation marginal effects plot
combined_precip_plot <- tmax_density / marg_precip +
  plot_layout(heights = c(2, 5), guides = "collect") +
  plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5))) &
  theme(plot.margin = margin(5, 5, 5, 5))


# Reduce gap between panels
marg_temp <- marg_temp +
  theme(plot.margin = margin(-5, 0, 5, 0))

marg_precip <- marg_precip +
  theme(plot.margin = margin(-5, 0, 5, 0))

# Combine panels
combined_temp_plot <- ppt_density / marg_temp +
  plot_layout(heights = c(1.5, 5)) & 
  theme(plot.margin = margin(0, 5, 0, 5))

combined_precip_plot <- tmax_density / marg_precip +
  plot_layout(heights = c(1.5, 5)) & 
  theme(plot.margin = margin(0, 5, 0, 5))

# Wrap plots for patchwork 
wrapped_temp_plot <- combined_temp_plot
wrapped_precip_plot <- combined_precip_plot

# Combine plots
final_figure <- wrapped_temp_plot | wrapped_precip_plot

final_figure
ggsave(here("Output", "fig5_marginal_effects.png"), 
       final_figure, 
       width = 12, height = 8, dpi = 300)

