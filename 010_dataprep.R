##########################################################################
### Author: Emma Nichols
### Project: LASIDAD Education - cog decline
### Purpose: Data Prep 
##########################################################################

rm(list = ls())

pacman::p_load(data.table, openxlsx, haven, readr, dplyr, magrittr, stringr,
               labelled, scales, gridExtra, grid, ggplot2, lubridate)
date <- gsub("-", "_", Sys.Date())
set.seed(6541)

# SET OBJECTS -------------------------------------------------------------

dropbox_dir <- "DIR"
dir <- paste0(dropbox_dir, "DIR")
lasi_raw_dir <- paste0(dropbox_dir, "DIR")
lasidad_w1_dir <- paste0(dropbox_dir, "DIR")
harmonized_dir <- paste0(dropbox_dir, "DIR")
longitudinal_dir <- paste0(dropbox_dir, "DIR")
biomarker_dir <- paste0(dropbox_dir, "DIR")
exit_dir <- paste0(dropbox_dir, "DIR")
rawdata_dir <- paste0(dir, "data/source/")
derived_dir <- paste0(dir, "data/derived/")

# GET DATA ----------------------------------------------------------------

lasi_data <- read_dta(paste0(harmonized_dir, "H_LASI_a3.dta"))
lasi_dt <- as.data.table(lasi_data)

state_map <- data.table(state_code = attr(lasi_data$hh1state, "labels"),
                        state_name = names(attr(lasi_data$hh1state, "labels")))
state_map[, state_name := gsub("[0-9]*\\.", "", state_name)]
state_map <- state_map[!is.na(state_code) & !state_code == 37] ## not missing or abroad

rawlasi_dt <- as.data.table(read_dta(paste0(lasi_raw_dir, "lasi_w1b_ind_bm.dta")))

lasidad_dt <- as.data.table(read_dta(paste0(harmonized_dir, "h_dad_w1a3.dta")))

exit_dt <- as.data.table(read_dta(paste0(exit_dir, "dad_exit_clean.dta"))) 

attrition_dt <- as.data.table(read_dta(paste0(rawdata_dir, "dad_tracker_v1.dta")))

longitudinal_dt <- as.data.table(read_dta(paste0(rawdata_dir, "H_LASI_DAD_b1.dta")))

varmap_dt <- as.data.table(read.xlsx(paste0(dir, "variable_map.xlsx")))

# FORMAT LONGITUDINAL DATA FOR SURVIVAL -----------------------------------------------------------

survival_dt <- copy(longitudinal_dt)
simple_vars <- varmap_dt[dataset == "longitudinal" & !grepl("\\*", variable), variable]
simple_vars_newname <- varmap_dt[dataset == "longitudinal" & variable %in% simple_vars, label]
long_vars <- varmap_dt[dataset == "longitudinal" & !variable %in% simple_vars, variable]

## RENAME SIMPLE VARS
setnames(survival_dt, simple_vars, simple_vars_newname)
survival_dt[, (simple_vars_newname) := lapply(.SD, as.numeric), .SDcols = simple_vars_newname] ## get rid of haven labelled

## RENAME LONG VARS
long_vars_old <- c(gsub("\\*", "1", long_vars), gsub("\\*", "2", long_vars))
long_vars_label <- varmap_dt[dataset == "longitudinal" & !variable %in% simple_vars, label]
long_vars_new <- c(paste0(long_vars_label, "1"), paste0(long_vars_label, "2"))
setnames(survival_dt, long_vars_old, long_vars_new)

## REMOVE REFRESHER SAMPLE
survival_dt <- survival_dt[!(inwave1 == 0 & inwave2 == 1)]

## CALCULATE FOLLOW-UP TIME FOR THOSE WHO CAME BACK
survival_dt[, date1 := my(paste0(iw_month1, "_", iw_year1))][, date2 := my(paste0(iw_month2, "_", iw_year2))]
survival_dt[, time_elapsed := as.period(date2-date1)/years(1)]

## ADD ON VARIABLES FOR W1 THAT ARE SIMPLE AND NEEDED FOR SURVIVAL DATA 
setnames(survival_dt, c("caste1", "age1", "rural1"), c("caste", "age", "rural"))

## SAVE VARIABLES NEEDED
survival_dt <- survival_dt[, c(simple_vars_newname, "caste", "age", "rural", "inwave2", "date1", "date2", "time_elapsed"), with = FALSE]

## CASTE CATEGORIZATIONS
survival_dt[, caste := factor(dplyr::case_when(caste == 1 ~ "Scheduled caste", 
                                      caste == 2 ~ "Scheduled tribe", 
                                      caste == 3 ~ "Other backward class", 
                                      caste == 4 ~ "No caste or other caste"), 
                     levels = c("No caste or other caste", "Scheduled caste", "Scheduled tribe", "Other backward class"))]

## AGE GROUPS
survival_dt[, age_group := as.numeric(cut(age, breaks = c(0,70,80,120), include.lowest = TRUE, right = FALSE))]
survival_dt[, age_group := factor(dplyr::case_when(age_group == 1 ~ "60-69", 
                                          age_group == 2 ~ "70-79", 
                                          age_group == 3 ~ "80+"), 
                         levels = c("60-69", "70-79", "80+"))]

survival_dt[, female := as.numeric(gender == 2)]                         

# CLASSIFIFY THOSE WHO REFUSED OR WERE SICK AS CENSORED AT VISIT DATE ------------

attrition_copy <- copy(attrition_dt)
attrition_copy <- attrition_copy[inw1 == 1]

attrition_copy[, `:=` (prim_key = as.numeric(prim_key), recontact_date = my(paste0(r2recontm, "_", r2reconty)), 
                       death_month = xd001_1, death_year = xd001_2,
                       died = r2died, refused = as.numeric(r2status %in% c(102,105)), attrited = as.numeric(r2status %in% c(103,104)))]

survival_dt <- merge(survival_dt, attrition_copy[, .(prim_key, recontact_date, death_month, death_year, died, attrited, refused)], all.x = TRUE)

## death date
int_midpoint <- function(interval){
    as.Date(interval@start + as.duration(interval)/2)
}
paste0("Only month missing: N = ", nrow(survival_dt[died == 1 & is.na(death_month) & !is.na(death_year)])) ## only month missing
survival_dt[died == 1 & is.na(death_month) & !is.na(death_year), death_month := 6] ## if month of death missing assume halfway through the year
survival_dt[, death_date := my(paste0(death_month, "_", death_year))]
paste0("Year or year and month missing: N = ", nrow(survival_dt[died == 1 & is.na(death_year)])) ## year of death missing
survival_dt[died == 1 & is.na(death_year), death_date := int_midpoint(interval(date1, recontact_date))] ## assume mid-point between W1 and recontact for those with missing year info
paste0("Incorrect (before W1 interview) death date: N = ", nrow(survival_dt[death_date < date1])) ## death date before W1 interview
survival_dt[died == 1 & death_date < date1, death_date := int_midpoint(interval(date1, recontact_date))] ## assume mid-point between W1 and recontact for those with incorrect death date

## date2 and time elapsed variables
survival_dt[refused == 1, date2 := recontact_date][died == 1, date2 := death_date]
survival_dt[refused == 1 | died == 1, time_elapsed := as.period(date2-date1)/years(1)]
paste0("Died in the same month as W1 interview: N = ", nrow(survival_dt[death_date == date1])) ## died in the same month as W1 interview
survival_dt[death_date == date1, time_elapsed := 1/24] ## add half a month for those who died in the same month as W1 interview

# FORMAT LONGITUDINAL DATA FOR LONGITUDINAL ANALYSIS --------------------------------

long_dt <- copy(longitudinal_dt)
simple_vars <- varmap_dt[dataset == "longitudinal" & !grepl("\\*", variable), variable]
simple_vars_newname <- varmap_dt[dataset == "longitudinal" & variable %in% simple_vars, label]
long_vars <- varmap_dt[dataset == "longitudinal" & !variable %in% simple_vars, variable]

## RENAME SIMPLE VARS
setnames(long_dt, simple_vars, simple_vars_newname)
long_dt[, (simple_vars_newname) := lapply(.SD, as.numeric), .SDcols = simple_vars_newname] ## get rid of haven labelled

## MAKE VARIABLES LONG
make_long <- function(var){
  wide_vars <- names(long_dt)[grepl(paste0(gsub("\\*", "[0-9]", var), "$"), names(long_dt))]
  long_chunk <- melt.data.table(long_dt, id.vars = c("prim_key"), measure.vars = wide_vars, 
                                variable.name = "wave", value.name = varmap_dt[variable == var, label])
  if (!var == "inw*"){ ## this works differently for the variable for being in a wave
    long_chunk[, wave := str_extract(wave, "^(h|r)[0-9]")][, wave := as.numeric(gsub("^.", "", wave))]
  }  else {
    long_chunk[, wave := str_extract(wave, "^inw[0-9]")][, wave := as.numeric(gsub("^inw", "", wave))]
  }
  if (long_chunk[, class(get(c(varmap_dt[variable == var, label])))][1] == "haven_labelled")
    long_chunk[, c(varmap_dt[variable == var, label]) := as.numeric(get(varmap_dt[variable == var, label]))] ## have to make numeric b/c haven labelled creates issues

  return(long_chunk)
}

long_chunks <- lapply(long_vars, make_long)
all_long <- Reduce(function(x,y) merge(x,y,all=TRUE,by=c("prim_key", "wave")), long_chunks)
long_dt <- merge(long_dt, all_long, by = "prim_key", all = TRUE)

## MARK REFRESHER SAMPLE
refresher_primkeys <- long_dt[wave == 1 & inwave == 0, unique(prim_key)]
long_dt[, refresher := as.numeric(prim_key %in% refresher_primkeys)]

## REMOVE THOSE NOT IN A WAVE
long_dt <- long_dt[!inwave == 0 & !is.na(inwave)]

## CREATE TIME VARIABLES
time_dt <- copy(survival_dt[died == 0 & attrited == 0 & refused == 0, .(prim_key, time = time_elapsed, wave = 2)])
long_dt <- merge(long_dt, time_dt, by = c("prim_key", "wave"), all = TRUE)
long_dt[wave == 1, time := 0]

## SAVE VARIABLES NEEDED
long_dt <- long_dt[, c(simple_vars_newname, varmap_dt[variable %in% long_vars, label], "time", "wave", "refresher"), with = FALSE]

## CASTE CATEGORIZATIONS
long_dt[, caste := factor(dplyr::case_when(caste == 1 ~ "Scheduled caste", 
                                      caste == 2 ~ "Scheduled tribe", 
                                      caste == 3 ~ "Other backward class", 
                                      caste == 4 ~ "No caste or other caste"), 
                     levels = c("No caste or other caste", "Scheduled caste", "Scheduled tribe", "Other backward class"))]

## AGE GROUPS
long_dt[, age_group := as.numeric(cut(age, breaks = c(0,70,80,120), include.lowest = TRUE, right = FALSE))]
long_dt[, age_group := factor(dplyr::case_when(age_group == 1 ~ "60-69", 
                                          age_group == 2 ~ "70-79", 
                                          age_group == 3 ~ "80+"), 
                         levels = c("60-69", "70-79", "80+"))]

## GENDER
long_dt[, female := as.numeric(gender == 2)]

## IW STATUS
long_dt[, status := factor(dplyr::case_when(status == 1 ~ "Cog + Inf", 
                                            status == 2 ~ "Cog only", 
                                            status == 3 ~ "Inf only"))]

## BINARY BY MEDIAN W1 COGNITION
w1_mediancog <- copy(long_dt[wave == 1])
median_cog <- w1_mediancog[, median(gcp, na.rm = TRUE)]
w1_mediancog[, above_w1medcog := as.numeric(gcp > median_cog)]
long_dt <- merge(long_dt, w1_mediancog[, .(prim_key, above_w1medcog)], by = "prim_key", all.x = TRUE)

## BINARY BY FOLLOW-UP TIME
w2_followup <- copy(long_dt[wave == 2])
w2_medfollowup <- w2_followup[, median(time, na.rm = TRUE)]
w2_followup[, above_w2medfu := as.numeric(time > w2_medfollowup)]
long_dt <- merge(long_dt, w2_followup[, .(prim_key, above_w2medfu)], by = "prim_key", all.x = TRUE)

## STANDARDIZE INDIVIDUAL TESTS BASED ON WAVE 1
tests <- c("wr_delayed", "lm_delayed", "animals", "con_praxis", "ravens")
for (var in tests){
    var_mean <- long_dt[wave == 1, mean(get(var), na.rm = TRUE)]
    var_sd <- long_dt[wave == 1, sd(get(var), na.rm = TRUE)]
    long_dt[, c(var) := (get(var)-var_mean)/var_sd]
}

## STANDARDIZE FUNCTION MEASURES TO WAVE 2 (FOR PRACTICE EFFECTS WEIGHTS)
ftests <- c("jorm", "blessed2", "blessed1")
for (var in ftests){
    var_mean <- long_dt[wave == 2, mean(get(var), na.rm = TRUE)]
    var_sd <- long_dt[wave == 2, sd(get(var), na.rm = TRUE)]
    long_dt[, c(var) := (get(var)-var_mean)/var_sd]
}

# PREP LASI VARIABLES ----------------------------------------------------------------

lasi_copy <- copy(lasi_dt)

## SET NAMES AND RESTRICT VARIABLES
setnames(lasi_copy, varmap_dt[dataset == "lasi", variable], varmap_dt[dataset == "lasi", label])
lasi_copy <- lasi_copy[, c("prim_key", varmap_dt[dataset == "lasi", label]), with = FALSE]

## CREATE WEALTH QUINTILE
lasi_copy[wealth < -quantile(wealth, probs = 0.05, na.rm = TRUE), wealth := NA] ## set to missing if more debt that the 5th percentile of wealth (likely have resources)
lasi_copy[, wealthq := as.numeric(cut(wealth, breaks = quantile(wealth, probs = seq(0,1,0.2), na.rm = TRUE), include.lowest = TRUE))]
lasi_copy[, wealthq := factor(case_when(wealthq == 1 ~ "Quintile 1",
                                        wealthq == 2 ~ "Quintile 2",
                                        wealthq == 3 ~ "Quintile 3",
                                        wealthq == 4 ~ "Quintile 4", 
                                        wealthq == 5 ~ "Quintile 5"))]

## CREATE CONSUMPTION QUINTILE
lasi_copy[, consumption := consumption_foodweek*52+consumption_nonfoodmonth*12+consumption_nonfoodyear]
lasi_copy[, consumptionq := as.numeric(cut(consumption, breaks = quantile(consumption, probs = seq(0,1,0.2), na.rm = TRUE), include.lowest = TRUE))]
lasi_copy[, consumptionq := factor(case_when(consumptionq == 1 ~ "Quintile 1",
                                             consumptionq == 2 ~ "Quintile 2",
                                             consumptionq == 3 ~ "Quintile 3",
                                             consumptionq == 4 ~ "Quintile 4", 
                                             consumptionq == 5 ~ "Quintile 5"))]
lasi_copy[, c("consumption_foodweek", "consumption_nonfoodmonth", "consumption_nonfoodyear") := NULL]

## EDUCATION VARIABLES
lasi_copy[, educ := factor(case_when(educ == 0 ~ "No school", 
                           educ == 1 ~ "Less than primary", 
                           educ == 2 ~ "Primary school", 
                           educ %in% 3:4 ~ "Middle-secondary school", 
                           educ > 4 ~ "Higher secondary school and up"), 
                    levels = c("No school", "Less than primary", "Primary school", 
                               "Middle-secondary school", "Higher secondary school and up"))]

peduc_vars <- c("educ_mom", "educ_dad")
lasi_copy[, (peduc_vars) := lapply(.SD, function(x) factor(case_when(x == 0 ~ "No school", 
                                                                    x %in% 1:2 ~ "Less than primary or primary", 
                                                                    x == 2 ~ "Primary school", 
                                                                    x > 3 ~ "Middle-secondary school and up"), 
                                                          levels = c("No school", "Less than primary or primary", "Middle-secondary school and up"))), .SDcols = peduc_vars]

## CHILDHOOD VARIABLES
lasi_copy[, childhood_finance := factor(case_when(childhood_finance == 1 ~ "Pretty well off", 
                                                  childhood_finance == 2 ~ "Average", 
                                                  childhood_finance == 3 ~ "Poor"), 
                                        levels = c("Average", "Poor", "Pretty well off"))]

lasi_copy[, childhood_health := factor(case_when(childhood_health == 1 ~ "Very good", 
                                                 childhood_health == 2 ~ "Good", 
                                                 childhood_health == 3 ~ "Fair", 
                                                 childhood_health == 4 ~ "Poor", 
                                                 childhood_health == 5 ~ "Very poor"), 
                                       levels = c("Fair", "Very poor", "Poor", "Good", "Very good"))]

## GET STATE NAME
lasi_copy <- merge(lasi_copy, state_map, by = "state_code")
lasi_copy[, statef := factor(state_name)]

## ADLS
lasi_copy[, adls := factor(case_when(adls == 0 ~ "None", 
                                     adls == 1 ~ "1", 
                                     adls == 2 ~ "2",
                                     adls %in% 3:4 ~ "3-4", 
                                     adls %in% 5:6 ~ "5-6"), 
                           levels = c("None", "1", "2", "3-4", "5-6"))]

## IADLS
lasi_copy[, iadls := factor(case_when(iadls == 0 ~ "None", 
                                      iadls == 1 ~ "1", 
                                      iadls == 2 ~ "2",
                                      iadls == 3 ~ "3",
                                      iadls %in% 4:5 ~ "4-5", 
                                      iadls %in% 6:7 ~ "6-7"), 
                            levels = c("None", "1", "2", "4-5", "6-7"))]

## STANDARDIZE GRIP STRENGTH 
lasi_copy[, gripstrength := (gripstrength-mean(gripstrength, na.rm = TRUE))/sd(gripstrength, na.rm = TRUE)]

## ELIGIABLE FOR W1
w1_ids <- long_dt[wave == 1, unique(prim_key)]
w1_states <- lasi_copy[prim_key %in% w1_ids, unique(state_code)]
lasi_copy[, elig_w1 := as.numeric(state_code %in% w1_states & age_lasi >=60)] ## eligible if in W1 states and age is greater than or equal to 60

lasi_copy[, prim_key := as.numeric(prim_key)] ## for merges

# PREP LASI-DAD W1 VARIABLES -----------------------------------------------------------

lasidad_copy <- copy(lasidad_dt)

## SET NAMES AND RESTRICT VARIABLES
setnames(lasidad_copy, varmap_dt[dataset == "lasidad_w1", variable], varmap_dt[dataset == "lasidad_w1", label])
lasidad_copy <- lasidad_copy[, c("prim_key", varmap_dt[dataset == "lasidad_w1", label]), with = FALSE]

## FORMAT CDR 
lasidad_copy[, cdr := factor(case_when(cdr == 0 ~ "Normal", 
                                       cdr == 0.5 ~ "MCI/Questionnable", 
                                       cdr >= 1 ~ "Dementia"), 
                             levels = c("Normal", "MCI/Questionnable", "Dementia"))]

lasidad_copy[, prim_key := as.numeric(prim_key)] ## for merges                 

# MERGE EVERYTHING TOGETHER -----------------------------------------------------------

survival_dt <- merge(survival_dt, lasi_copy, by = "prim_key", all.x = TRUE)
survival_dt <- merge(survival_dt, lasidad_copy, by = "prim_key", all.x = TRUE)

long_dt <- merge(long_dt, lasi_copy, by = "prim_key", all.x = TRUE)
long_dt <- merge(long_dt, lasidad_copy, by = "prim_key", all.x = TRUE)

# MAKE PRACTICE EFFECTS DATA ---------------------------------------------------------

pe_dt <- copy(long_dt)
nrow(pe_dt[wave == 2 & elig_w1 == 0])
pe_dt <- pe_dt[wave == 2 & elig_w1 == 1]
nrow(pe_dt[status == "Inf only"])
pe_dt <- pe_dt[!status == "Inf only"] ## remove those with only informant data

# REMOVE REFRESHER SAMPLE FOR LONGITUDINAL DATA --------------------------------------

long_dt <- long_dt[refresher == 0]

# CLEAN UP DATA ----------------------------------------------------------------------

## GET RID OF HAVEN LABELED
for (var in names(survival_dt)){
    if (survival_dt[, class(get(var))][1] == "haven_labelled") survival_dt[, c(var) := as.numeric(get(var))]
}
for (var in names(long_dt)){
    if (long_dt[, class(get(var))][1] == "haven_labelled") long_dt[, c(var) := as.numeric(get(var))]
}
for (var in names(pe_dt)){
    if (pe_dt[, class(get(var))][1] == "haven_labelled") pe_dt[, c(var) := as.numeric(get(var))]
}

write_rds(list(survival = survival_dt, longitudinal = long_dt, practice = pe_dt), paste0(derived_dir, "processed_data.rds"))
