##==============================================================================
# Purpose: Create SDTM DS dataset for Genentech ADS Programmer assignment
#
# Auth: DLW
# Date: 7/14/2026
#
# Reference documents
#   - eCRF: https://github.com/pharmaverse/pharmaverseraw/blob/main/vignettes/articles/aCRFs/Subject_Disposition_aCRF.pdf
#.  - CT: https://raw.githubusercontent.com/pharmaverse/examples/refs/heads/main/metadata/sdtm_ct.csv
#   - Raw DS: pharmaverseraw::ds_raw
#.  - Raw DM: pharmaverseraw::dm_raw
#
# Notes: 
# An error-free program with good documentation that will create the DS domain 
# with the following variables: STUDYID, DOMAIN, USUBJID, DSSEQ, DSTERM, DSDECOD, 
# DSCAT, VISITNUM, VISIT, DSDTC, DSSTDTC, DSSTDY
# 
# Deliverables
# - SDTM creation script: question_1_sdtm/01_create_ds_domain.R
# - Resulting SDTM dataset in any format
# - A text file/log file as evidence for code running error-free
##==============================================================================

library(dplyr) 
library(sdtm.oak)
library(admiral)
library(lubridate) # remove 
library(tidyr) # remove?

## Get necessary data
# load DS raw data
ds_raw <- pharmaverseraw::ds_raw

# load DM raw data (for RFICDTC var)
dm_raw <- pharmaverseraw::dm_raw

# read in control terminology
ct.url <- "https://raw.githubusercontent.com/pharmaverse/examples/refs/heads/main/metadata/sdtm_ct.csv"
ct <- read.csv(ct.url, stringsAsFactors = F)

## 
# all blanks look to already be NA
ds_raw <- convert_blanks_to_na(dat)

# create oak id var
ds_raw <- ds_raw |>
  generate_oak_id_vars(pat_var = "PATNUM",
                       raw_src = "ds_raw")

## data checks
# IT.DSDECOD and OTHERSP should never both be complete
table(is.na(ds_raw$OTHERSP), is.na(ds_raw$IT.DSDECOD))
# IT.DSDECOD and OTHERSP use ct C66727
ct$collected_value

# generate ds domain
ds <- 
  # derive topic variable
  # map IT.DSTERM to DSTERM if OTHERSP is NA, otherwise map OTHERSP to DSTERM
  # DSTERM = if_else(is.na(OTHERSP), IT.DSTERM, OTHERSP)
  assign_no_ct(
    raw_dat = condition_add(ds_raw, !is.na(ds_raw$OTHERSP)),
    raw_var = "OTHERSP",
    tgt_var = "DSTERM",
    id_vars = oak_id_vars()
  ) |>
  assign_no_ct(
    raw_dat = condition_add(ds_raw, is.na(ds_raw$OTHERSP)),
    raw_var = "IT.DSTERM",
    tgt_var = "DSTERM",
    id_vars = oak_id_vars()
  ) |>
  # DSDECOD has controlled terminology
  # map OTHERSP to DSDECOD if OTHERSP is non-missing, otherwise map IT.DSDECOD to DSDECOD
  assign_ct(
    raw_dat = condition_add(ds_raw, !is.na(OTHERSP)),
    raw_var = "OTHERSP", 
    tgt_var = "DSDECOD",
    ct_spec = ct,
    ct_clst = "C66727",
    id_vars = oak_id_vars()
  ) |>
  assign_ct(
    raw_dat = condition_add(ds_raw, is.na(OTHERSP)),
    raw_var = "IT.DSDECOD", 
    tgt_var = "DSDECOD",
    ct_spec = ct,
    ct_clst = "C66727",
    id_vars = oak_id_vars()
  ) |>
  # DSCAT 
  hardcode_ct(
    raw_dat = condition_add(ds_raw, IT.DSDECOD == "Randomized"),
    raw_var = "IT.DSDECOD",
    tgt_var = "DSCAT",
    tgt_val = "PROTOCOL MILESTONE",
    ct_clst = "C74558",
    ct_spec = ct,
  ) |>
  hardcode_ct(
    raw_dat = condition_add(ds_raw, IT.DSDECOD != "Randomized"),
    raw_var = "IT.DSDECOD",
    tgt_var = "DSCAT",
    tgt_val = "DISPOSITION EVENT",
    ct_clst = "C74558",
    ct_spec = ct,
  ) |>
  hardcode_ct(
    raw_dat = condition_add(ds_raw, !is.na(OTHERSP)),
    raw_var = "OTHERSP",
    tgt_var = "DSCAT",
    tgt_val = "OTHER EVENT",
    ct_clst = "C74558",
    ct_spec = ct
  ) |>
  # DSSTDTC
  assign_datetime(
    raw_dat = ds_raw, 
    raw_var = "IT.DSSTDAT",
    tgt_var = "DSSTDTC",
    raw_fmt = "m-d-y"
  ) |>
  # DSDTC
  assign_datetime(
    raw_dat = ds_raw,
    raw_var = c("DSDTCOL", "DSTMCOL"),
    tgt_var = "DSDTC",
    raw_fmt = c("m-d-y", "H:M")
  ) |>
  mutate(DOMAIN = "DS") |>
  # Visit information
  assign_ct(
    raw_dat = ds_raw,
    raw_var = "INSTANCE",
    tgt_var = "VISIT",
    ct_spec = ct,
    ct_clst = "VISIT"
  ) |>
 assign_ct(
   raw_dat = ds_raw,
   raw_var = "INSTANCE",
   tgt_var = "VISITNUM",
   ct_spec = ct,
   ct_clst = "VISITNUM"
 )

|>
  # DSSEQ
  derive_seq(
    rec_vars = c("VISIT", "DSTERM"),
    sbj_vars = c("STUDY", "PATNUM"),
    tgt_var = "DSSEQ"
  )


# DPLYR version 
# need to add mappings to controlled terminology
ds_dplyr <- ds_raw |>
  rename(STUDYID = STUDY,
         USUBJID = PATNUM) |>
  mutate(DOMAIN = "DS",
         DSDECOD = if_else(is.na(OTHERSP), IT.DSDECOD, OTHERSP),
         DSCAT = case_when(IT.DSDECOD == "Randomized" ~ "PROTOCOL MILESTONE",
                           !is.na(IT.DSDECOD) ~ "DISPOSITION EVENT",
                           !is.na(OTHERSP) ~ "OTHER EVENT"),
         DSTERM = if_else(is.na(OTHERSP), IT.DSTERM, OTHERSP),
         DSSTDTC = format_ISO8601(as.Date(ds.dplyr$IT.DSSTDAT, format = "%m-%d-%Y")))
        # need to add another date variable





