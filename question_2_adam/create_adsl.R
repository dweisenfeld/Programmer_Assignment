##==============================================================================
# Purpose: Create mock ADaM ADSL dataset for Genentech ADS Programmer assignment
#
# Auth: DLW
# Date: 7/14/2026
#
# Input data
#   - pharmaversesdtm::dm
#   - pharmaversesdtm::vs
#   - pharmaversesdtm::ex
#   - pharmaversesdtm::ds 
#   - pharmaversesdtm::ae
#
# Output:
#   - ADSL dataset saved as RDS: "output/adsl.rds"
#
# Notes: 
# - derive these variables
#   - AGEGR9/AGEGR9N:   <19, 18-50, >50 / Age group number (1,2,3)
#   - TRTSDTM/TRTSTMF:  Earliest datetime w/ valid dose / time imputation flag
#   - ITTFL:            intention-to-treat flag (Y/N)
#   - LSTAVLDT:         Last date known to be alive 
##==============================================================================

sink("logs/question2_log.txt")

library(admiral)
library(dplyr)
library(labelled)

out_path <- "./output/"

# load data
dm <- pharmaversesdtm::dm
vs <- pharmaversesdtm::vs
ex <- pharmaversesdtm::ex
ds <- pharmaversesdtm::ds 
ae <- pharmaversesdtm::ae

################################# Pre-Process ################################## 
# standardize blanks (this isn't necessary here, but is for SAS data)
dm <- convert_blanks_to_na(dm)
vs <- convert_blanks_to_na(vs)
ex <- convert_blanks_to_na(ex)
ds <- convert_blanks_to_na(ds)
ae <- convert_blanks_to_na(ae)

# EX - add time component to date field
ex_ext <- ex |>
 derive_vars_dtm(
   dtc              = EXSTDTC,
   new_vars_prefix  = "EXST",
   time_imputation  = "00:00:00",
   flag_imputation  = "time"
 ) |>
 derive_vars_dtm(
   dtc              = EXENDTC,
   new_vars_prefix  = "EXEN",
   time_imputation  = "00:00:00",
   flag_imputation  = "time"
 ) |> 
 # create flag to identify valid doses
 mutate(valid_dose_flag = (EXDOSE > 0 | (EXDOSE==0 & EXTRT=="PLACEBO"))) 

################################# Generate ADSL ################################

adsl <- dm |> select(-DOMAIN) 

## Create lookups for categorical variables
# Age group
agegr9_lookup <- exprs(
  ~condition,               ~AGEGR9,     ~AGEGR9N,
  AGE < 18,                 "<18",       1,
  AGE >= 18 & AGE <= 50,    "18 - 50",   2,
  AGE > 50,                 ">50",       3
)
# ITT Flag
ittfl_lookup <- exprs(
  ~condition,     ~ITTFL,
  !is.na(ARM),    "Y",
  is.na(ARM),     "N"
)

adsl <- adsl |>
  # age group
  derive_vars_cat(
    definition = agegr9_lookup
  ) |>
  # ITT flag
  derive_vars_cat(
   definition = ittfl_lookup
  ) |>
  # 
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = valid_dose_flag==T, # should condition be here rather than creating a flag?
    new_vars = exprs(TRTSDTM = EXSTDTM, TRTSTMF = EXSTTMF),
    order = exprs(EXSTDTM, EXSEQ),
    mode = "first",
    by_vars = exprs(STUDYID, USUBJID)
  ) |>
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = valid_dose_flag==T, 
    new_vars = exprs(TRTEDTM = EXENDTM, TRTETMF = EXENTMF),
    order = exprs(EXENDTM, EXSEQ),
    mode = "last",
    by_vars = exprs(STUDYID, USUBJID)
  ) |>
  derive_vars_dtm_to_dt(source_vars = exprs(TRTSDTM, TRTEDTM))

# LSTAVLDT
# Set to the last date patient has documented clinical data to show him/her
# alive, converted to numeric date, using the following dates:
# (1) last complete date of vital assessment with a valid test result
#      ([VS.VSSTRESN] and [VS.VSSTRESC] not both missing) and datepart of
#      [VS.VSDTC] not missing.
# (2) last complete onset date of AEs (datepart of Start Date/Time of
#     Adverse Event [AE.AESTDTC]).
# (3) last complete disposition date (datepart of Start Date/Time of
#     Disposition Event [DS.DSSTDTC]).
# (4) last date of treatment administration where patient received a valid
#     dose (datepart of Datetime of Last Exposure to Treatment
#     [ADSL.TRTEDTM]).
# Set to max of (Vitals complete, AE onset complete, disposition complete,
#                treatment complete).

# Use derive_vars_extreme_event to derive LSTAVLDT
adsl <- adsl |> 
  derive_vars_extreme_event(
    by_vars = exprs(STUDYID, USUBJID),
    events = list(
      # last complete date of vital assessment with a valid test result
      event(
        dataset_name = "vs",
        order = exprs(VSDTC, VSSEQ),
        condition = !is.na(VSSTRESN) | !is.na(VSSTRESC),
        set_values_to = exprs(
          LSTAVLDT = convert_dtc_to_dt(VSDTC, highest_imputation = "M"),
          seq = VSSEQ
        ),
      ),
      # last complete onset date of AEs
      event(
        dataset_name = "ae",
        order = exprs(AESTDTC, AESEQ),
        condition = !is.na(AESTDTC),
        set_values_to = exprs(
          LSTAVLDT = convert_dtc_to_dt(AESTDTC, highest_imputation = "M"),
          seq = AESEQ
        ),
      ),
      # last complete disposition date
      event(
        dataset_name = "ds",
        order = exprs(DSSTDTC, DSSEQ),
        condition = !is.na(DSSTDTC),
        set_values_to = exprs(
          LSTAVLDT = convert_dtc_to_dt(DSSTDTC, highest_imputation = "M"),
          seq = DSSEQ
        ),
      ),
      # last date of treatment administration where patient received a valid dose
      event(
        dataset_name = "adsl",
        condition = !is.na(TRTEDT),
        set_values_to = exprs(LSTAVLDT = TRTEDT, seq = 0),
      )
    ),
    source_datasets = list(ae = ae, vs = vs, ds = ds, adsl = adsl),
    tmp_event_nr_var = event_nr,
    order = exprs(LSTAVLDT, seq, event_nr),
    mode = "last",
    new_vars = exprs(LSTAVLDT)
  )

adsl.final <- adsl |>
  select(
    STUDYID, USUBJID, SUBJID, SITEID,
    AGE, AGEU, AGEGR9N, AGEGR9,
    SEX, RACE, ETHNIC, ARM, ACTARM, ITTFL,
    TRTSDTM, TRTSTMF, TRTEDTM, TRTETMF, TRTEDT,
    LSTAVLDT,
    everything()
  ) |>
  arrange(USUBJID)

## Add labels for derived variables
var_label(adsl.final) <- list(
  AGEGR9     = "Age Group 9",
  AGEGR9N    = "Age Group 9 (N)",
  TRTSDTM    = "Datetime of First Exposure to Treatment",
  TRTSTMF    = "Time of First Exposure Imputation Flag",
  ITTFL      = "Intent-To-Treat Population Flag",
  LSTAVLDT   = "Last Date Known Alive"
)


############################### Data checks ####################################
# is there one record per subject
nrow(adsl.final)==length(unique(adsl$USUBJID))
nrow(adsl.final)==length(unique(dm$USUBJID))

# check age mappings
table(adsl$AGEGR9, adsl$AGEGR9N, exclude = NULL) 

# check variable missingness
adsl |>
  summarise(across(c(AGEGR9, AGEGR9N, TRTSDTM, TRTSTMF, ITTFL, LSTAVLDT), ~ sum(is.na(.x)))) 
  # only derived TRT vars have missingness and those should for screen failures


################################## Save ########################################
# Save the derived ADSL dataset as RDS
# could use xportr to write a CDISC-compliant .xpt file if metadata was available
saveRDS(adsl.final, file = paste0(out_path, "adsl.rds"))

sink()
