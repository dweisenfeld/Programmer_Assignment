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

################################# Load data ####################################
# load DS raw data
ds_raw <- pharmaverseraw::ds_raw

# load DM raw data (for RFICDTC var)
dm_raw <- pharmaverseraw::dm_raw
  # IC_DT maps to RFICDTC

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

############################## data checks #####################################
get_terms <- function(code) {
  unique(unlist(ct[ct$codelist_code==code, c("collected_value", "term_preferred_term", "term_synonyms")]))
}

# IT.DSDECOD and OTHERSP should never both be complete
table(is.na(ds_raw$OTHERSP), is.na(ds_raw$IT.DSDECOD))

# IT.DSDECOD and OTHERSP use ct C66727
dsdecod_terms <- get_terms("C66727")
setdiff(c(ds_raw$OTHERSP, ds_raw$IT.DSDECOD), dsdecod_terms)
setdiff(dsdecod_terms, c(ds_raw$OTHERSP, ds_raw$IT.DSDECOD))
  # there are minor inconsistencies in many of the raw values that don't
  # properly map to CT

# CT for VISIT
visit_terms <- get.terms("VISIT")
setdiff(ds_raw$INSTANCE, visit_terms)
setdiff(visit_terms, ds_raw$INSTANCE)

# CT for VISITNUM
visitnum_terms <- get.terms("VISITNUM")
setdiff(ds_raw$INSTANCE, visitnum_terms)
setdiff(visitnum_terms, ds_raw$INSTANCE)
intersect(visitnum_terms, ds_raw$INSTANCE)
  # most VISIT/VISITNUM options map to CT, except for unscheduled visits. 
  # Only Unscheduled 3.1 maps, need to map other unscheduled visits


# adding additional rows to CT to account for minor inconsistencies
ct_updates <- tribble(
 ~codelist_code, ~term_code,  ~term_value,              ~collected_value,
 "VISIT",        "VISIT",     "Ambul ECG Removal",      "Ambul Ecg Removal",
 "VISIT",        "VISIT",     "UNSCHEDULED 1.1",        "Unscheduled 1.1",
 "VISIT",        "VISIT",     "UNSCHEDULED 4.1",        "Unscheduled 4.1",
 "VISIT",        "VISIT",     "UNSCHEDULED 5.1",        "Unscheduled 5.1",
 "VISIT",        "VISIT",     "UNSCHEDULED 6.1",        "Unscheduled 6.1",
 "VISIT",        "VISIT",     "UNSCHEDULED 8.2",        "Unscheduled 8.2",
 "VISIT",        "VISIT",     "UNSCHEDULED 13.1",       "Unscheduled 13.1",
 "VISITNUM",     "VISITNUM",  "Ambul ECG Removal",      "Ambul Ecg Removal",
 "VISITNUM",     "VISITNUM",  "UNSCHEDULED 1.1",        "Unscheduled 1.1",
 "VISITNUM",     "VISITNUM",  "UNSCHEDULED 4.1",        "Unscheduled 4.1",
 "VISITNUM",     "VISITNUM",  "UNSCHEDULED 5.1",        "Unscheduled 5.1",
 "VISITNUM",     "VISITNUM",  "UNSCHEDULED 6.1",        "Unscheduled 6.1",
 "VISITNUM",     "VISITNUM",  "UNSCHEDULED 8.2",        "Unscheduled 8.2",
 "VISITNUM",     "VISITNUM",  "UNSCHEDULED 13.1",       "Unscheduled 13.1" 
)
  
# adding term synonyms to correct minor differences between actual collected value
# and values listed in ct
update
ct_updated <- ct |> bind_rows(ct_updates)
ct_updated$term_synonyms[ct_updated$collected_value=="Lost To Follow-Up"] <- "Lost to Follow-Up"
ct_updated$term_synonyms[ct_updated$collected_value=="Study Terminated By Sponsor"] <- "Study Terminated by Sponsor"
ct_updated$term_synonyms[ct_updated$collected_value=="Trial Screen Failure"] <- "Screen Failure"
ct_updated$term_synonyms[ct_updated$collected_value=="Complete"] <- "Completed"


############################## generate ds domain ##############################
ds <- 
  # create STUDYID, USUBJID, and DOMAIN
  assign_no_ct(
    raw_dat = ds_raw,
    raw_var = "STUDY",
    tgt_var = "STUDYID"
  ) |>
  assign_no_ct(
    raw_dat = ds_raw,
    raw_var = "PATNUM",
    tgt_var = "USUBJID"
  ) |>
  mutate(DOMAIN = "DS") |>
  # derive topic variable
  # map IT.DSTERM to DSTERM if OTHERSP is NA, otherwise map OTHERSP to DSTERM
  # DSTERM = if_else(is.na(OTHERSP), IT.DSTERM, OTHERSP)
  assign_no_ct(
    raw_dat = condition_add(ds_raw, !is.na(ds_raw$OTHERSP)),
    raw_var = "OTHERSP",
    tgt_var = "DSTERM"
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
    ct_spec = ct_updated,
    ct_clst = "C66727",
    id_vars = oak_id_vars()
  ) |>
  assign_ct(
    raw_dat = condition_add(ds_raw, is.na(OTHERSP)),
    raw_var = "IT.DSDECOD", 
    tgt_var = "DSDECOD",
    ct_spec = ct_updated,
    ct_clst = "C66727",
    id_vars = oak_id_vars()
  ) |>
  # DSCAT 
  # if IT.DSDECOD = Randomized, DSCAT = PROTOCOL MILESTONE, else DSCAT = DISPOSITION EVENT
  # if OTHERSP is non-missing then DSCAT = OTHER EVENT
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
  # Visit information (has CT)
  assign_ct(
    raw_dat = ds_raw,
    raw_var = "INSTANCE",
    tgt_var = "VISIT",
    ct_spec = ct_updated,
    ct_clst = "VISIT"
  ) |>
 assign_ct(
   raw_dat = ds_raw,
   raw_var = "INSTANCE",
   tgt_var = "VISITNUM",
   ct_spec = ct_updated,
   ct_clst = "VISITNUM"
 ) |>
  # DSSEQ
  derive_seq(
    rec_vars = c("VISIT", "VISITNUM"),
    sbj_vars = c("STUDYID", "USUBJID"), 
    tgt_var = "DSSEQ"
  ) |>
  select(
    STUDYID, 
    DOMAIN, 
    USUBJID, 
    DSSEQ, 
    DSTERM, 
    DSDECOD, 
    DSCAT, 
    VISITNUM, 
    VISIT, 
    DSDTC, 
    DSSTDTC
  )

## derive study day
# modify DM domain raw data 
dm <- dm_raw |>
  mutate(RFSTDTC = create_iso8601(IC_DT, .format = "m/d/y"),
         USUBJID = PATNUM)

ds <- ds |>
  derive_study_day(
    dm_domain = dm, 
    tgdt = "DSSTDTC",  
    refdt = "RFSTDTC", 
    study_day_var = "DSSTDY"
)

# there are still 3 terms that could not be mapped per CT
#     - OTHERSP: Final Lab Visit, Final Retrieval Visit
#     - IT.DSDECOD: Randomized


### DPLYR version 
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





