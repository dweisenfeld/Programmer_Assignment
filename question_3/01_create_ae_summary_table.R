##==============================================================================
# Purpose: Create and save summary table of treatment-emergent AEs (TEAEs)
#
# Auth: DLW
# Date: 7/19/2026
#
# Input datasets
#   - pharmaverseadam::adae 
#.  - pharmaverseadam::adsl
# 
# Deliverables
# - Script to create summary table: question_3_tlg/01_create_ae_summary_table.R
# - Text files/log files as evidence for code running error-free
# - ae_summary_table.html (or .docx/.pdf)
#
# Notes:
# Summary Table using {gtsummary} - HINT - FDA Table 10
# Create a summary table of treatment-emergent adverse events (TEAEs).
# - Treatment-emergent AE records will have TRTEMFL == "Y" in pharmaverseadam::adae
# - Rows: AETERM or AESOC
# - Columns: Treatment groups (ACTARM)
# - Cell values: Count (n) and percentage (%)
# - Include total column with all subjects
# - Sort by descending frequency
##==============================================================================

sink("logs/question3.01_log.txt")

library(gtsummary)
library(gt)
library(dplyr)

out_path <- "./output/"

## Load data
adae <- pharmaverseadam::adae 
adsl <- pharmaverseadam::adsl

## Pre-process data
teae <- adae |>
  filter(TRTEMFL == "Y") |> # treatment emergent AE records 
  select(USUBJID, AETERM, AESOC, ACTARM)

## Create table  
teae_tbl <- teae |>
  tbl_hierarchical(
    variables = c(AESOC, AETERM),
    by = ACTARM,
    id = USUBJID,
    denominator = adsl,
    overall_row = TRUE,
    label = "..ard_hierarchical_overall.." ~ "Any SAE"
  ) |>
  sort_hierarchical(sort = "descending")

## Save table
teae_tbl |> 
  as_gt() |> 
  gtsave(filename = paste0(out_path, "/ae_summary_table.html"))

sink()
