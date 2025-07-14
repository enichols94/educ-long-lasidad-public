##########################################################################
### Author: Emma Nichols
### Date: 01/29/2025
### Project: LASIDAD Educationa and longitudinal change
### Purpose: RTM simulation
##########################################################################

rm(list = ls())

pacman::p_load(data.table, openxlsx, haven, readr, dplyr, magrittr, stringr,
               labelled, scales, gridExtra, grid, ggplot2, lubridate, JM, 
               survival, patchwork, survey, srvyr, doBy, lme4, pbapply, future.apply)
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
plot_dir <- paste0(dir, "paper/mortality_fig/")

iterations <- 1000

# READ DATA ------------------------------------------------------------------

data <- read_rds(paste0(derived_dir, "processed_data_weights.rds"))
survival_dt <- data$survival
longitudinal_dt <- data$longitudinal

# RESTRICTIONS ---------------------------------------------------------------

## number lost to mortality and loss to follow-up  
survival_dt[, sum(died)]; survival_dt[, sum(as.numeric(refused == 1 | attrited == 1))]

## remove missing data
dt <- copy(longitudinal_dt); surv_dt <- copy(survival_dt)
model_covs <- c("age", "gender", "caste", "rural", "educ_dad", "childhood_finance", "childhood_health")
for (var in c("gcp", "educ", model_covs)){
    message(paste0("Missing individuals in ", var, ": ", dt[is.na(get(var)), length(unique(prim_key))], " (", 
                   dt[is.na(get(var)) & wave == 2, length(unique(prim_key))], " longitudinal)"))
    dt <- dt[!is.na(get(var))]
    message(paste0("Now excluding N=", nrow(surv_dt[!prim_key %in% dt[, unique(prim_key)]]), " from data completely"))
}
surv_dt <- surv_dt[prim_key %in% dt[, unique(prim_key)]] ## only keep those in longitudinal data

# ADJUST DATA ----------------------------------------------------------------

## any school binary variable
survival_dt[, any_school := as.numeric(!educ == "No school")]
longitudinal_dt[, any_school := as.numeric(!educ == "No school")]

## add gcp to survival data
survival_dt <- merge(survival_dt, longitudinal_dt[wave == 1, .(prim_key, gcp)], by = "prim_key", all.x = TRUE)

# SIMULATION -----------------------------------------------------------------

## exposure: any school binary
## w1 cognition: based on education
## w2 cognition: defined by w1 cognition * any school

# MODELS TO INFORM SIM --------------------------------------------------------

## prevalence of any school 
anyschool_mean <- survival_dt[, mean(any_school)]

## w1 cog
w1cog_model <- lm(gcp ~ any_school, data = longitudinal_dt[wave == 1])
w1cog_params <- as.data.table(parameters::model_parameters(w1cog_model))

## w2 cog
w2cog_dt <- na.omit(dcast.data.table(longitudinal_dt, prim_key + any_school ~ wave, value.var = "gcp", fill = NA))
setnames(w2cog_dt, c("1", "2"), c("w1_cog", "w2_cog"))
w2cog_model <- lm(w2_cog ~ w1_cog, data = w2cog_dt)
w2cog_params <- as.data.table(parameters::model_parameters(w2cog_model))

# CREATE SIM FUNCTION -----------------------------------------------------------

sim_data <- function(sample_size = 2000, w2_effect){

    ## create data
    sim_dt <- data.table(prim_key = 1:sample_size)
    sim_dt[, any_school := rbinom(sample_size, size = 1, prob = anyschool_mean)]
    sim_dt[, w1_gcp := w1cog_params[Parameter == "(Intercept)", Coefficient] + 
              w1cog_params[Parameter == "any_school", Coefficient] * any_school + 
              rnorm(sample_size, 0, sigma(w1cog_model))] # add noise to w1 cognition    
    sim_dt[, w2_gcp := w2cog_params[Parameter == "(Intercept)", Coefficient] + 
              w2_effect * w1_gcp + 
              rnorm(sample_size, 0, sigma(w2cog_model))]

    ## results to check sim parameterization 
    prop_anyschool <- sim_dt[, mean(any_school)]
    cog_anyschool <- parameters::model_parameters(lm(w1_gcp ~ any_school, data = sim_dt))$Coefficient[2]
    w2gcp_w1gcp <- parameters::model_parameters(lm(w2_gcp ~ w1_gcp, data = sim_dt))$Coefficient[2]
    
    ## reformat longitudinal
    lsim_dt <- melt.data.table(sim_dt, id.vars = c("prim_key", "any_school"), measure.vars = c("w1_gcp", "w2_gcp"), 
                               variable.name = "wave", value.name = "gcp")
    lsim_dt[, wave := as.numeric(str_extract(wave, "[0-9]"))]
    lsim_dt[, time := ifelse(wave == 1, 0, 5)]

    ## run models with and without mortality
    model_sim <- lmer(gcp ~ time * any_school + (1 | prim_key), data = lsim_dt)
    params_sim <- as.data.table(parameters::model_parameters(model_sim))

    ## format results
    params_select <- grepl("time:any_school", params_sim$Parameter)
    result_dt <- data.table(w2w1_param = w2_effect,
                            int_value = params_sim$Coefficient[params_select], 
                            prop_anyschool = prop_anyschool, cog_anyschool = cog_anyschool, w2gcp_w1gcp = w2gcp_w1gcp)
    return(result_dt)
}

full_sim_results <- data.table()
for (effect in seq(0.4,1.6,0.2)){
    message(paste0("Running simulation with w2 effect: ", effect))
    sim_results <- rbindlist(pbreplicate(iterations, 
                                sim_data(w2_effect = effect), 
                                simplify = FALSE))
    full_sim_results <- rbind(full_sim_results, sim_results)
}

# PRINT SIMULATION RESULTS -----------------------------------------------------------

## check simulation parameterization
full_sim_results[, lapply(.SD, mean), by = w2w1_param, .SDcols = c("prop_anyschool", "cog_anyschool", "w2gcp_w1gcp")] 

# PLOT RESULTS ------------------------------------------------------------------------

sim_summary <- full_sim_results[, .(effect = w2w1_param, mean = mean(int_value), lower = quantile(int_value, 0.025), upper = quantile(int_value, 0.975)), 
                                by = "w2w1_param"]

sim_plot <- ggplot(sim_summary, aes(x = as.factor(effect), y = mean, ymin = lower, ymax = upper)) + 
    geom_point() +
    geom_errorbar(width = 0) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(x = "Coefficient estimate (effect of W1 on W2)", y = "Difference in cognitive slope (any education vs. none)") +
    scale_color_manual(name = "", values = paletteer::paletteer_c("grDevices::Temps", 7)) + 
    theme_bw()

write_rds(sim_plot, paste0(derived_dir, "simulated_rtm_plot.rds"))
