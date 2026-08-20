# Run simulation 2 for the regression-to-the-mean sensitivity analysis.
#
# Input: data/derived/processed_data_weights.rds from 011_weights.R.
# Output: data/derived/simulation2.rds for script 052.

source(file.path("R", "project_paths.R"))
load_required_packages(c(
  "data.table", "openxlsx", "haven", "readr", "dplyr", "magrittr",
  "stringr", "labelled", "scales", "gridExtra", "grid", "ggplot2",
  "lubridate", "JM", "paletteer", "survival", "patchwork", "survey",
  "srvyr", "doBy", "lme4", "pbapply", "future.apply", "parameters"
))
set.seed(6541)

# PATHS AND SIMULATION SETTINGS ------------------------------------------

paths <- project_paths()
derived_dir <- paths$derived
ensure_output_directory(derived_dir)

iterations <- 1000

# READ DATA ------------------------------------------------------------------

data <- read_rds(required_input_file(file.path(derived_dir, "processed_data_weights.rds")))
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

## w2 true cog
change_dt <- dcast.data.table(longitudinal_dt, prim_key ~ wave, value.var = "gcp", fill = NA)
change_dt <- change_dt[!is.na(`2`)]
change_dt[, change := `2` - `1`]
mean_change <- change_dt[, mean(change)]
sd_change <- change_dt[, sd(change)]

# CREATE SIM FUNCTION -----------------------------------------------------------

sim_data <- function(sample_size = 2000, prop.noise = 0.8){

    ## create data
    sim_dt <- data.table(prim_key = 1:sample_size)
    sim_dt[, any_school := rbinom(sample_size, size = 1, prob = anyschool_mean)]
    sim_dt[, w1_true := w1cog_params[Parameter == "(Intercept)", Coefficient] + 
              w1cog_params[Parameter == "any_school", Coefficient] * any_school + 
              rnorm(sample_size, 0, sigma(w1cog_model)*(1-prop.noise))]
    sim_dt[, w1_obs := w1_true + rnorm(sample_size, 0, sigma(w1cog_model)*prop.noise)] 
    sim_dt[, true_change := rnorm(sample_size, mean = mean_change, sd = sd_change*(1-prop.noise))]
    sim_dt[, w2_obs := w1_true + true_change + ## if you change w1_true to w1_obs then you get no correlation
              rnorm(sample_size, 0, sd_change*prop.noise)] 
    sim_dt[, w2_true := w1_true + true_change]      
    sim_dt[, obs_change := w2_obs - w1_obs]    

    ## reformat longitudinal - observed
    lsim_dt <- melt.data.table(sim_dt, id.vars = c("prim_key", "any_school"), measure.vars = c("w1_obs", "w2_obs"), 
                               variable.name = "wave", value.name = "gcp")
    lsim_dt[, wave := as.numeric(str_extract(wave, "[0-9]"))]
    lsim_dt[, time := ifelse(wave == 1, 0, 5)]

    ## reformat longitudinal - true
    lsim_true <- melt.data.table(sim_dt, id.vars = c("prim_key", "any_school"), measure.vars = c("w1_true", "w2_true"), 
                                  variable.name = "wave", value.name = "gcp")
    lsim_true[, wave := as.numeric(str_extract(wave, "[0-9]"))]
    lsim_true[, time := ifelse(wave == 1, 0, 5)]

    ## run models
    #model_obs <- lmer(gcp ~ time * any_school + (1 | prim_key), data = lsim_dt)
    model_obs <- lm(gcp ~ time * any_school, data = lsim_dt) ## will give same mean answer as linear model because of data generation 
    params_obs <- as.data.table(parameters::model_parameters(model_obs))

    ## get correlations
    obs_cor <- sim_dt[, cor(w1_obs, obs_change)]

    ## format results
    params_select <- grepl("time:any_school", params_obs$Parameter)
    result_dt <- data.table(noise_proportion = prop.noise, type = c("observed"),
                            int_value = params_obs$Coefficient[params_select], 
                            change_school = sim_dt[any_school == 1, mean(obs_change)],
                            change_noschool = sim_dt[any_school == 0, mean(obs_change)],
                            correlation = c(obs_cor))
    return(result_dt)
}

full_sim_results <- data.table()
for (noise in c(1)){ ## doesn't really matter what the noise proportion is, so only show 1
    message(paste0("Running simulation with noise proportion: ", noise))
    sim_results <- rbindlist(pbreplicate(iterations, 
                                sim_data(prop.noise = noise), 
                                simplify = FALSE))
    full_sim_results <- rbind(full_sim_results, sim_results)
}

# PRINT SIMULATION RESULTS -----------------------------------------------------------


# PLOT RESULTS ------------------------------------------------------------------------

sim_summary <- full_sim_results[, .(mean = mean(int_value), lower = quantile(int_value, 0.025), upper = quantile(int_value, 0.975)), 
                                by = "noise_proportion"]
sim_summary[, noise_proportion := "Estimate"]                                

sim2_plot <- ggplot(sim_summary, aes(x = as.factor(noise_proportion), y = mean, ymin = lower, ymax = upper)) + 
    geom_point() +
    geom_errorbar(width = 0) + 
    lims(y = c(-0.08, 0.08)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(x = "", y = "Difference in cognitive slope \n(any education vs. none)") +
    scale_color_manual(name = "", values = paletteer::paletteer_c("grDevices::Temps", 7)) + 
    theme_bw()

write_rds(sim2_plot, file.path(derived_dir, "simulation2.rds"))
