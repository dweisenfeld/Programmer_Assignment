# Programmer Assignment

ADS programmer assignment demonstrating data programming for clinical trials in R using the `sdtm.oak` and `admiral` packages. 

## Overview of each script

### Question 1 — SDTM DS Domain (`question_1_sdtm/`)
Maps raw disposition data (`pharmaverseraw::ds_raw`) to a SDTM DS domain using `sdtm.oak` and saves the output to `output/ds.rds`.  

### Question 2 — ADaM ADSL (`question_2_adam/`)
Derives a mock ADaM `ADSL` dataset from SDTM sources (data from `pharmaversesdtm`) using `admiral` and saves the output to `output/adsl.rds`.

### Question 3 — AE Summary & Visualizations (`question_3/`)
Creates tables and figures to summarize study AEs
- **`01_create_ae_summary_table.R`** — builds a treatment-emergent adverse event (TEAE) summary table by treatment arm with `gtsummary`, saved as `output/ae_summary_table.html`.
- **`02_create_visualizations.R`** — generates a bargraph showing AE severity by arm (`output/AEs_by_Severity_Arm.png`) and a figure showing the top 10 AEs by Organ Class (`output/Top_AEs.png`).

## Requirements

This project uses R and RStudio. Required packages:

```r
install.packages(c("dplyr", "labelled", "gtsummary", "gt", "ggplot2", "binom", "scales"))

# pharmaverse packages
install.packages(c("admiral", "sdtm.oak", "sdtmchecks"))
install.packages(c("pharmaverseraw", "pharmaversesdtm", "pharmaverseadam"))
```

## Running the Scripts

Open `project.Rproj` in RStudio, then run the scripts from the repository root. Scripts can be run in any order. 

Each script writes a console log to `logs/` and saves its resulting dataset, table, or plot to `output/`.

## Output

Script | File | Description |
|---|---|---|
`question_1_sdtm/X` | `output/ds.rds` | SDTM DS domain |
`question_2_adam/X` | `output/adsl.rds` | ADaM ADSL dataset |
`question_3/01_create_ae_summary_table.R` | `output/ae_summary_table.html` | TEAE summary table |
`question_3/02_create_visualizations.R` | `output/AEs_by_Severity_Arm.png` | AE severity by treatment arm |
`question_3/02_create_visualizations.R` | `output/Top_AEs.png` | Top 10 AEs with 95% CIs |
