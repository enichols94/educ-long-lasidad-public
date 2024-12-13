##########################################################################
### Author: Emma Nichols
### Date: 12/11/2024
### Project: LASIDAD Educationa and longitudinal change
### Purpose: Descriptive analyses
##########################################################################

rm(list = ls())

pacman::p_load(data.table, openxlsx, haven, readr, dplyr, magrittr, stringr,
               labelled, scales, gridExtra, grid, ggplot2, lubridate, JM, 
               survival, patchwork, survey, srvyr)
date <- gsub("-", "_", Sys.Date())
set.seed(6541)

# SET OBJECTS -------------------------------------------------------------

dropbox_dir <- "C:/Users/emmanich/P2AGING Dropbox/Emma Nichols/"
dir <- paste0(dropbox_dir, "projects/educ_long_lasidad/")
lasi_raw_dir <- paste0(dropbox_dir, "H_LASI/ToUpload/Raw/Data/LASI_w1b_Stata/")
harmonized_dir <- paste0(dropbox_dir, "Harmonized Data Files/")
longitudinal_dir <- paste0(dropbox_dir, "H_DAD/Raw_wave2/Preliminary LASI-DAD-Core/")
exit_dir <- paste0(dropbox_dir, "H_DAD/Raw_wave2/Combined/Data/Clean/")
rawdata_dir <- paste0(dir, "data/source/")
derived_dir <- paste0(dir, "data/derived/")
plot_dir <- paste0(dir, "plots/")

# READ DATA ------------------------------------------------------------------

data <- read_rds(paste0(derived_dir, "processed_data_weights.rds"))
survival_dt <- data$survival
longitudinal_dt <- data$longitudinal

dt <- copy(survival_dt)
dt <- merge(dt, longitudinal_dt[wave == 1, .(prim_key, gcp, memory, executive, language, visuospatial)], by = "prim_key")

design <- dt[!died == 1 & !refused == 1 & !attrited == 1] %>% as_survey_design(ids = 1, weights = full_weight)

# CREATE OUTCOME VARIABLE -----------------------------------------------------

dt[, outcome := factor(case_when(inwave2 == 1 ~ "In wave 2", 
                                 refused == 1 ~ "Refused wave 2", 
                                 died == 1 ~ "Died between waves", 
                                 attrited == 1 ~ "Lost to follow-up"), 
                       levels = c("In wave 2", "Refused wave 2", "Died between waves", "Lost to follow-up"))]

# CREATE TABLE ----------------------------------------------------------------

## weighted/unweighted with combo weight (sampling, death, attrition)

## outcome, follow-up time, age, gender, caste, rural, educ_mom, educ_dad, childhood_finance_childhood_health
## improvedsanitation, puccahouse, wealthq, consumptionq

get_col <- function(educ_cat){

    if (educ_cat == "Overall"){
        d <- design; data <- copy(dt)
    } else {
        d <- design %>% subset(educ == educ_cat)
        data <- copy(dt[educ == educ_cat])  
    }

    ## HELPER FUNCTIONS
    get_continuous <- function(char, num = 1){
        west <- d %>% srvyr::summarize(mean = survey_mean(get(char), na.rm = TRUE))
        est <- data[, mean(get(char), na.rm = TRUE)]
        wsd <- d %>% srvyr::summarize(sd = survey_sd(get(char), na.rm = TRUE))
        sd <- data[, sd(get(char), na.rm = TRUE)]
        paste0(sprintf(paste0("%.", num, "f"), est), " (", sprintf(paste0("%.", num, "f"), sd), ") / ",
               sprintf(paste0("%.", num, "f"), west[1]), " (", sprintf(paste0("%.", num, "f"), wsd[1]), ")")
    }

    get_binary <- function(char, cat){
        west <- d %>% srvyr::summarize(mean = survey_mean(get(char)==cat, proportion = T, na.rm = T))*100
        est <- nrow(data[get(char) == cat])/nrow(data)*100
        N <- nrow(data[get(char) == cat])
        paste0(sprintf("%.1f", est), " / ", sprintf("%.1f", west[1]), " (", nrow(data[get(char)==cat]), ")")
    } 

    get_binary_raw <- function(char, cat){
        paste0(sprintf("%.1f", nrow(data[get(char) == cat])/nrow(data)*100), " (", nrow(data[get(char) == cat]), ")")
    }

    get_continuous_raw <- function(char, num = 1, df = data){
        paste0(sprintf(paste0("%.", num, "f"), df[, mean(get(char), na.rm = TRUE)]), " (", 
               sprintf(paste0("%.", num, "f"), df[, sd(get(char), na.rm = TRUE)]), ")")
    }

    ## GET VALUES
    outcome <- sapply(dt[, levels(outcome)], function(x) get_binary_raw("outcome", x))
    followup <- get_continuous_raw("time_elapsed", df = data[inwave2 == 1])
    
    age <- get_continuous("age")
    female <- get_binary("female", 1)
    rural <- get_binary("rural", 1)
    caste <- sapply(dt[, levels(caste)], function(x) get_binary("caste", x))

    educ_mom <- sapply(dt[, levels(educ_mom)], function(x) get_binary("educ_mom", x))
    educ_dad <- sapply(dt[, levels(educ_dad)], function(x) get_binary("educ_dad", x))
    child_finance <- sapply(dt[, levels(childhood_finance)], function(x) get_binary("childhood_finance", x))
    child_health <- sapply(dt[, levels(childhood_health)], function(x) get_binary("childhood_health", x))

    col <- data.table(name = c("", outcome, followup, age, female, rural, "", caste, "", educ_mom, "", educ_dad,
                               "", child_finance, "", child_health))
    setnames(col, "name", educ_cat)
    return(col)
}

all_cols <- lapply(c("Overall", dt[, levels(educ)]), get_col)
table_dt <- Reduce(function(x,y) cbind(x,y), all_cols)

labels <- c("Wave 2 status", paste0("  ", dt[, levels(outcome)]), "Follow-up time", "Age (yrs)", "Female (%)", "Rural (%)", "Caste", 
            paste0("  ", dt[, levels(caste)]), "Maternal education", paste0("  ", dt[, levels(educ_mom)]), "Paternal education", 
            paste0("  ", dt[, levels(educ_dad)]), "Childhood financial status", paste0("  ", dt[, levels(childhood_finance)]), 
            "Childhood health status", paste0("  ", dt[, levels(childhood_health)]))

table_dt <- cbind(labels, table_dt)
write.xlsx(table_dt, paste0(plot_dir, "table1_", date, ".xlsx"))
