##==============================================================================
# Purpose: Create and save summary table of treatment-emergent AEs (TEAEs)
#
# Auth: DLW
# Date: 7/19/2026
#
# Input data
#   - pharmaverseadam::adae 
#   - pharmaverseadam::adsl 
#
# Outputs:
# - Plot 1: AEs_by_Severity_Arm.png
# - Plot 2: Top_AEs.png
#
# Notes:
# - Plot 1: AE severity distribution by treatment (bar chart or heatmap). 
#   AE Severity is captured in the AESEV variable in pharmaverseadam::adae dataset.
# - Plot 2: Top 10 most frequent AEs (with 95% CI for incidence rates). 
#   AEs are captured inthe AETERM variable in the pharmaverseadam::adae dataset.
##==============================================================================

sink("logs/question3.02_log.txt")

library(dplyr)
library(ggplot2)
library(binom)     # for Clopper-Pearson CIs
library(scales)    

fig_path <- "./output/"

## Load data
adae <- pharmaverseadam::adae 
adsl <- pharmaverseadam::adsl

##################### PLOT 1 - AE SEVERITY BY TREATMENT ########################

aesev_counts <- adae |> 
  group_by(ARM, AESEV) |>
  summarise(count = n(), .groups = "drop_last") 

aesev_counts |>
  ggplot(aes(x = ARM, y = count, fill = AESEV)) +
  geom_col()

ae_sev_barplot <- adae |>
  mutate(y = 1) |> # create variable for ggplot to sum up
  ggplot(aes(x = ARM, y = y, fill = AESEV)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "bottom") +
  labs(y = "Count of AEs", x = "Treatment Arm",
       title = "Adverse Event Severity Distribution by Treatment Arm",
       fill = "AE Severity")

# save figure
ggsave(filename = paste0(fig_path, "AEs_by_Severity_Arm.png"), 
       plot = ae_sev_barplot,
       width = 6, height = 5, units = "in")

##### PLOT 2 - Top 10 most frequent AEs (with 95% CI for incidence rates) ######

# need to determine total num of subjects in safety population 
n_subjects <- adsl |>
  filter(SAFFL=="Y") |>
  n_distinct("USUBJID") 

# calculate count of each AE (only counted once per subject)
ae_counts <- adae |>
  distinct(USUBJID, AETERM) |>    # one record per subject per AE
  summarise(ae_count = n(), .by = AETERM)

# Calculate Clopper-Pearson 95% 
ae_ci <- binom.confint(
  x = ae_counts$ae_count,
  n = n_subjects, # or length(unique(adae$USUBJID)) to make example?
  conf.level = 0.95,
  methods = "exact" 
) 

# create plot data -- need AETERM in order of frequency
ae_plot_df <- ae_ci |>
  bind_cols(ae_counts) |>
  arrange(-mean) |>
  slice_head(n = 10) |>
  mutate(AETERM = factor(AETERM, levels = rev(AETERM)))  

# plot
top_aes_plot <- ae_plot_df |>
  ggplot(aes(x = mean, y = AETERM)) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = .4) +
  geom_point(size = 4) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Top 10 Most Frequent Adverse Events",
       subtitle = paste0("n = ", n_subjects, " subjects; 95% Clopper-Pearson CIs"),
       x = "Percentage of Patients (%)",
       y = "") +
  theme_classic(base_size = 12) 

# save figure
ggsave(filename = paste0(fig_path, "Top_AEs.png"), top_aes_plot,
       width = 8, height = 5, units = "in")

## Note:
# The example figure shows N = 225 subjects, which is the number of subjects 
# with any AE. I based the denominator for calculating the incident rates on the
# total number of subjects in the safety population, which is 254, not 225.
# If it was instead desired to only use the number of subjects with AEs, then
# in the binom.confint function n_subjects should be replaced with 
# length(unique(adae$USUBJID))

sink(type = "message")
