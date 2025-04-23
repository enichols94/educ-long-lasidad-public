##########################################################################
### Author: Emma Nichols
### Date: 01/17/2025
### Project: LASIDAD Educationa and longitudinal change
### Purpose: Practice effects analysis
##########################################################################

rm(list = ls())

pacman::p_load(data.table, openxlsx, haven, readr, dplyr, magrittr, stringr,
               labelled, scales, gridExtra, grid, ggplot2, lubridate, JM, 
               survival, patchwork, survey, srvyr, doBy)
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
plot_dir <- paste0(dir, "paper/practiceeffects_fig/")

# READ DATA ------------------------------------------------------------------

pe_dt <- read_rds(paste0(derived_dir, "processed_practiceeffects_data.rds"))

cog_dt <- as.data.table(read_dta(paste0(rawdata_dir, "H_DAD_w12_noimp_core.dta")))

pe_map <- as.data.table(read.xlsx(paste0(dir, "practiceeffects_variable_map.xlsx")))

# FORMATTING FOR PRACTICE EFFECTS DATA ---------------------------------------

practice_dt <- copy(cog_dt)

simple_vars <- pe_map[!grepl("\\*", variable), variable]
simple_vars_newname <- pe_map[variable %in% simple_vars, label]
long_vars <- pe_map[!variable %in% simple_vars, variable]

## FIX INCONSISTENT VARIABLES - MAKE CONSISTENT ACROSS WAVES (despite administration differences)
setnames(practice_dt, c("r1sc_scorep", "r2sc_scoret", "r1hammer1", "r2hammer2"), 
                      c("r1sc_score", "r2sc_score", "r1hammer", "r2hammer"))

## RENAME SIMPLE VARS
setnames(practice_dt, simple_vars, simple_vars_newname)
practice_dt[, (simple_vars_newname) := lapply(.SD, as.numeric), .SDcols = simple_vars_newname] ## get rid of haven labelled

## MAKE VARIABLES LONG
make_long <- function(var){
  print(var)
  wide_vars <- names(practice_dt)[grepl(paste0(gsub("\\*", "[0-9]", var), "$"), names(practice_dt))]
  long_chunk <- melt.data.table(practice_dt, id.vars = c("prim_key"), measure.vars = wide_vars, 
                                variable.name = "wave", value.name = pe_map[variable == var, label])
  if (!var == "inw*"){ ## this works differently for the variable for being in a wave
    long_chunk[, wave := str_extract(wave, "^(h|r)[0-9]")][, wave := as.numeric(gsub("^.", "", wave))]
  }  else {
    long_chunk[, wave := str_extract(wave, "^inw[0-9]")][, wave := as.numeric(gsub("^inw", "", wave))]
  }
  if (long_chunk[, class(get(c(pe_map[variable == var, label])))][1] == "haven_labelled")
    long_chunk[, c(pe_map[variable == var, label]) := as.numeric(get(pe_map[variable == var, label]))] ## have to make numeric b/c haven labelled creates issues

  return(long_chunk)
}

long_chunks <- lapply(long_vars, make_long)
all_long <- Reduce(function(x,y) merge(x,y,all=TRUE,by=c("prim_key", "wave")), long_chunks)
practice_dt <- merge(practice_dt[, ..simple_vars_newname], all_long, by = "prim_key", all = TRUE)

## MERGE WAVE 2 ONTO PRACTICE EFFECTS DATA 
practice_dt <- practice_dt[wave == 2 & status %in% 1:2]
practice_dt[, c("status", "wave") := NULL]
pe_dt[, c("wr_delayed", "lm_delayed", "animals", "con_praxis", "ravens") := NULL] ## get rid of all cognitive tests (with imputations)
practice_dt <- merge(pe_dt, practice_dt, by = "prim_key", all.x = TRUE)

## CREATE NAMING VARIABLE 
practice_dt[, naming := scissors + elbow + hammer + pencil + watch]

## CREATE MISSING INDICATOR FOR LASI COGNITION
practice_dt[, gcp_lasiM := as.numeric(is.na(gcp_lasi))]
practice_dt[gcp_lasiM == 1, gcp_lasi := 0]

## BINARY EDUCATION VARIABLE 
practice_dt[, any_school := as.numeric(!educ == "No school")]

## MAKE SURE ALL VARIABLES ARE Z-SCORED
cog_vars <- c("gcp", "orientation", "memory", "executive", "language", "visuospatial", 
              "time_orient", "place_orient", "serial7", "naming", "wr_immediate", "wr_delayed", "lm_delayed", "animals", 
              "con_praxis", "ravens", "symbolcancellation", "gonogo", "tokentest")
for (var in cog_vars){
    practice_dt[, paste0(var, "Z") := scale(get(var), center = TRUE, scale = TRUE)]
}

# ESTIMATE PRACTICE EFFECTS -----------------------------------------------------

covs <- c("age", "female", "caste", "rural", "gcp_lasi", "gcp_lasiM", "educ")
cog_varsZ <- paste0(cog_vars, "Z")

get_practiceeffect <- function(var){

    ## copy data
    data <- copy(practice_dt)
    data_m <- copy(data)
    variables <- c(var, "practice", "pe_weight", covs) 

    ## format data
    data <- na.omit(data[, ..variables])
    data_m[, paste0(var, "_m") := as.numeric(is.na(get(var)))]
    data_design <- svydesign(ids = ~1, data = data, weights = ~pe_weight)
    data_m_design <- svydesign(ids = ~1, data = data_m, weights = ~pe_weight)

    ## model
    formula <- paste0(var, " ~ practice + ", paste(covs, collapse = " + "))
    model <- svyglm(as.formula(formula), design = data_design)
    params <- parameters::model_parameters(model)
    select_params <- grepl("practice", params$Parameter)
    results <- data.table(var = gsub("Z", "", var), 
                          estimate = params$Coefficient[select_params], 
                          lower = params$CI_low[select_params], 
                          upper = params$CI_high[select_params], 
                          sample_size = nrow(data), 
                          missingness_test = 0)

    ## missingness model if exists
    if (nrow(data_m[get(paste0(var, "_m")) == 1]) > 0){
        formulaM <- paste0(var, "_m ~ practice + ", paste(covs, collapse = " + "))
        modelM <- svyglm(as.formula(formulaM), design = data_m_design, family = binomial(link = "logit"))
        paramsM <- parameters::model_parameters(modelM)
        select_paramsM <- grepl("practice", paramsM$Parameter)
        resultsM <- data.table(var = gsub("Z", "", var), 
                               estimate = paramsM$Coefficient[select_paramsM], 
                               lower = paramsM$CI_low[select_paramsM], 
                               upper = paramsM$CI_high[select_paramsM], 
                               sample_size = nrow(data_m), 
                               missingness_test = 1)
    } else {
        resultsM <- data.table(var = gsub("Z", "", var), 
                               estimate = NA, 
                               lower = NA, 
                               upper = NA, 
                               sample_size = NA, 
                               missingness_test = 1)
    }

    results <- rbind(results, resultsM)
    return(results)
    
}

result_dt <- rbindlist(lapply(cog_varsZ, get_practiceeffect))

# ESTIMATE INTERACTIONS -----------------------------------------------------------

get_practiceeffect_strat <- function(var){

    ## copy data
    data <- copy(practice_dt)
    variables <- c(var, "practice", "pe_weight", "any_school", covs[!covs %in% "educ"]) 

    ## format data
    data <- na.omit(data[, ..variables])
    data_design <- svydesign(ids = ~1, data = data, weights = ~pe_weight)

    ## model
    formula <- paste0(var, " ~ practice*any_school + ", paste(covs[!covs %in% "educ"], collapse = " + "))
    model <- svyglm(as.formula(formula), design = data_design)
    params <- parameters::model_parameters(model)

    ## get stratified results
    coef_nums <- matrix(c(as.numeric(grepl("^practice$", params$Parameter)), 
                          as.numeric(grepl("practice", params$Parameter))),
                        byrow = TRUE, nrow = 2)
    combs <- doBy::esticon(model, L = coef_nums, conf.int = TRUE)

    select_int <- grepl("practice:any_school", params$Parameter)
    results <- data.table(var = gsub("Z", "", var), educ = c("No school", "Any school"),
                          estimate = combs$estimate,  
                          lower = combs$lwr, 
                          upper = combs$upr, 
                          p_int = params$p[select_int],
                          sample_size = nrow(data), 
                          missingness_test = 0)

    return(results)
    
}

results_strat_dt <- rbindlist(lapply(cog_varsZ, get_practiceeffect_strat))

# GRAPH RESULTS ------------------------------------------------------------------

varlabel_dt <- data.table(var = c("gcp", "orientation", "memory", "executive", "language", "visuospatial", 
                                  "time_orient", "place_orient", "serial7", "naming", "wr_immediate", "wr_delayed", "lm_delayed", 
                                  "animals", "con_praxis", "ravens", "symbolcancellation", "gonogo", "tokentest"),
                          label = c("Global cognitive performance", "Orientation", "Memory", "Executive function", "Language", "Visuospatial", 
                                    "Time orientation", "Place orientation", "Serial 7s", "Naming", "Word recall immediate", "Word recall delayed", "Logical memory delayed", 
                                    "Animals", "Constructional praxis", "Raven's matrices", "Symbol cancellation", "Go/no-go", "Token test"))
varlabel_dt[, label := factor(label, levels = rev(varlabel_dt[, label]))]                                    

## main data
plot_dt <- copy(result_dt)
plot_dt <- merge(plot_dt, varlabel_dt, by = "var")

## interaction data
plot_strat_dt <- copy(results_strat_dt)
plot_strat_dt <- merge(plot_strat_dt, varlabel_dt, by = "var")
plot_strat_dt[, educ := factor(educ, levels = rev(c("No school", "Any school")))]

## results section 
plot_dt[var == "gcp" & missingness_test == 0]
plot_dt[var %in% c("memory", "visuospatial", "lm_delayed") & missingness_test == 0]

## main results
main_graph <- ggplot(plot_dt[missingness_test == 0], aes(x = label, y = estimate, ymin = lower, ymax = upper)) + 
    geom_point(color = "#2972b6") + 
    geom_errorbar(width = 0.1, color = "#2972b6") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() + 
    labs(x = "", y = "Practice effect\n(Difference in standardized score between those\nwith practice vs. without)") + 
    theme_bw() 

## missingness results
breaks <- c(0.5,0.75,1,1.25,1.5,1.75,2)
log_breaks <- log(breaks)
missingness_graph <- ggplot(plot_dt[missingness_test == 1], aes(x = label, y = estimate, ymin = lower, ymax = upper)) + 
    geom_point(color = "#2972b6") + 
    geom_errorbar(width = 0.1, color = "#2972b6") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_y_continuous(breaks = log_breaks, labels = breaks, limits = c(-0.5, 0.7), oob = oob_keep) +
    coord_flip() + 
    labs(x = "", y = "Practice effect\n(Odds ratio for missing data comparing those\nwith practice vs. without)") + 
    theme_bw() 

## interaction results
educ_graph <- ggplot(plot_strat_dt, aes(x = label, y = estimate, ymin = lower, ymax = upper, color = educ)) + 
    geom_point(position = position_dodge(width = 0.5)) + 
    geom_errorbar(width = 0.1, position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(name = "Education", values = c("#8ec582", "#008083")) +
    coord_flip() + 
    guides(color = guide_legend(reverse = TRUE)) +
    labs(x = "", y = "Practice effect\n(Difference in standardized score between those\nwith practice vs. without)") + 
    theme_bw() + 
    theme(legend.position = "bottom")

## full graph 
full_graph <- (main_graph + missingness_graph) / 
              (plot_spacer() + educ_graph + plot_spacer() + plot_layout(widths = c(0.3, 2, 1.5))) +
              plot_layout(nrow = 2) +
              plot_annotation(tag_levels = "A", tag_suffix = ".")

ggsave(paste0(plot_dir, "practice_effects_", date, ".pdf"), plot = full_graph, width = 12, height = 8)
