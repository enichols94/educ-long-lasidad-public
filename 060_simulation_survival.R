##########################################################################
### Author: Emma Nichols
### Project: LASIDAD Education and longitudinal change
### Purpose: Mortality simulation
##########################################################################

rm(list = ls())

pacman::p_load(data.table, openxlsx, haven, readr, dplyr, magrittr, stringr,
               labelled, scales, gridExtra, grid, ggplot2, lubridate, JM, 
               survival, patchwork, survey, srvyr, doBy, lme4, pbapply, future.apply)
date <- gsub("-", "_", Sys.Date())
set.seed(6541)

# SET OBJECTS -------------------------------------------------------------

dropbox_dir <- "DIR"
dir <- paste0(dropbox_dir, "DIR")
lasi_raw_dir <- paste0(dropbox_dir, "DIR")
harmonized_dir <- paste0(dropbox_dir, "DIR")
longitudinal_dir <- paste0(dropbox_dir, "DIR")
exit_dir <- paste0(dropbox_dir, "DIR")
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
surv_dt[, any_school := as.numeric(!educ == "No school")]
dt[, any_school := as.numeric(!educ == "No school")]

## add gcp to survival data
surv_dt <- merge(surv_dt, dt[wave == 1, .(prim_key, gcp)], by = "prim_key", all.x = TRUE)

# SIMULATION -----------------------------------------------------------------

## exposure: any school binary
## w1 cognition: based on education
## w2 cognition: defined by w1 cognition * any school
## mortality: defined by w1 cognition * any school

# MODELS TO INFORM SIM --------------------------------------------------------

## prevalence of any school 
anyschool_mean <- surv_dt[, mean(any_school)]

## survival
surv_model <- glm(died ~ gcp*any_school, data = surv_dt, family = binomial(link = "logit"))
survmodel_params <- as.data.table(parameters::model_parameters(surv_model))

## w1 cog
w1cog_model <- lm(gcp ~ any_school, data = dt[wave == 1])
w1cog_params <- as.data.table(parameters::model_parameters(w1cog_model))

## w2 cog
w2cog_dt <- na.omit(dcast.data.table(dt, prim_key + any_school ~ wave, value.var = "gcp", fill = NA))
setnames(w2cog_dt, c("1", "2"), c("w1_cog", "w2_cog"))
w2cog_model <- lm(w2_cog ~ w1_cog, data = w2cog_dt)
w2cog_params <- as.data.table(parameters::model_parameters(w2cog_model))

# PLOT PROBABILITY OF DEATH -----------------------------------------------------

death_effects <- as.data.table(effects::effect("gcp*any_school", surv_model, xlevels = list(gcp=seq(-2,2,0.5), any_school = c(0,1))))

death_extra <- surv_dt[, .(gcp = mean(gcp)), by = "any_school"]
death_extra[, death_prob := predict.glm(surv_model, newdata = death_extra, type = "response")]
death_extra

death_plot <- ggplot() +
    geom_line(data = death_effects, aes(x = gcp, y = fit, color = as.factor(any_school))) + 
    geom_ribbon(data = death_effects, aes(x = gcp, y = fit, ymin = lower, ymax = upper, fill = as.factor(any_school)), alpha = 0.1, color = NA) +
    geom_vline(data = death_extra, aes(xintercept = gcp, color = as.factor(any_school)), linetype = "dashed") +
    geom_point(data = death_extra, aes(x = gcp, y = death_prob, color = as.factor(any_school)), size = 2) +
    labs(x = "Cognition", y = "Probability of death") +
    scale_color_manual(name = "Education", values = c("#ffb93c", "#ff585d"), labels = c("None", "Any")) +
    scale_fill_manual(name = "Education", values = c("#ffb93c", "#ff585d"), labels = c("None", "Any")) +
    theme_bw() + 
    theme(legend.position = "none")

deathplot_bottom <- ggplot(surv_dt, aes(x = gcp, fill = as.factor(any_school))) +
    geom_density(alpha = 0.5) +
    labs(x = "Cognition", y = "Density") +
    scale_x_continuous(limits = c(-2,2)) +
    scale_fill_manual(name = "Education", values = c("#ffb93c", "#ff585d"), labels = c("None", "Any")) +
    theme_bw() + 
    theme(legend.position = "bottom")

deathplot_full <- death_plot + deathplot_bottom + 
    plot_layout(ncol = 1, heights = c(4, 1))

ggsave(paste0(plot_dir, "death_plot_", date, ".pdf"), plot = deathplot_full, width = 8, height = 6)

# CALIBRATE INTERCEPTS ----------------------------------------------------------

## function to optimize over 
interceptopt <- function(intercept, gcp_effect, anyschool_effect, interaction_effect){
    
    data <- copy(surv_dt)

    set.seed(12345)

    data[, prop_died := plogis(intercept + 
                               gcp_effect * gcp + 
                               anyschool_effect * any_school + 
                               interaction_effect * gcp * any_school)]
    data[, sim_died := rbinom(nrow(data), size = 1, prob = prop_died)]
  
    optSreturn<-abs(data[, sum(sim_died)] - data[, sum(died)])
    
    return(optSreturn)
}

## do optimization
optimize_intercept <-function(nrep=500, Lwin=-3, Uwin=-1.5, g, s, i){
    
  set.seed(12345)
  
  plan(multisession, workers=10) ## Run in parallel on local computer. Can specify 10 cores for desktop.
  
  beta_S0_reps<-future_replicate(nrep,optimize(interceptopt, lower = Lwin, upper = Uwin,
                                               maximum = FALSE,
                                               gcp_effect = g, anyschool_effect = s, interaction_effect = i)$minimum)

  plan(sequential)
  
  return(mean(beta_S0_reps))
} 

intercept_main <- optimize_intercept(nrep=500, Lwin = -2, Uwin = -1, g = survmodel_params[Parameter == "gcp", Coefficient], 
                                     s = survmodel_params[Parameter == "any_school", Coefficient], 
                                     i = survmodel_params[Parameter == "gcp:any_school", Coefficient])

intercept_extreme <- optimize_intercept(nrep=500, Lwin = -2, Uwin = 0, g = log(0.5), 
                                         s = log(0.5), i = log(0.5))


# CREATE SIM FUNCTION -----------------------------------------------------------

sim_data <- function(sample_size = 2000, intercept = intercept_main, 
                     gcp_effect = survmodel_params[Parameter == "gcp", Coefficient],
                     anyschool_effect = survmodel_params[Parameter == "any_school", Coefficient],
                     interaction_effect = survmodel_params[Parameter == "gcp:any_school", Coefficient]){

    ## create data
    sim_dt <- data.table(prim_key = 1:sample_size)
    sim_dt[, any_school := rbinom(sample_size, size = 1, prob = anyschool_mean)]
    sim_dt[, w1_gcp := w1cog_params[Parameter == "(Intercept)", Coefficient] + 
              w1cog_params[Parameter == "any_school", Coefficient] * any_school + 
              rnorm(sample_size, 0, sigma(w1cog_model))]
    sim_dt[, w2_gcp := w2cog_params[Parameter == "(Intercept)", Coefficient] + 
              w2cog_params[Parameter == "w1_cog", Coefficient] * w1_gcp + 
              rnorm(sample_size, 0, sigma(w2cog_model))]
    sim_dt[, prop_died := plogis(intercept + 
                                 gcp_effect * w1_gcp + 
                                 anyschool_effect * any_school + 
                                 interaction_effect * w1_gcp * any_school)]
    sim_dt[, died := rbinom(sample_size, size = 1, prob = prop_died)]

    ## results to check sim parameterization 
    prop_anyschool <- sim_dt[, mean(any_school)]
    cog_anyschool <- parameters::model_parameters(lm(w1_gcp ~ any_school, data = sim_dt))$Coefficient[2]
    w2gcp_w1gcp <- parameters::model_parameters(lm(w2_gcp ~ w1_gcp, data = sim_dt))$Coefficient[2]
    surv_params <- parameters::model_parameters(glm(died ~ w1_gcp * any_school, data = sim_dt, family = binomial(link = "logit")))
    died_w1gcp <- surv_params$Coefficient[2]; died_anyschool <- surv_params$Coefficient[3]; died_int <- surv_params$Coefficient[4]
    num_died <- sim_dt[died == 1, .N]

    ## reformat longitudinal
    lsim_dt <- melt.data.table(sim_dt, id.vars = c("prim_key", "any_school", "died"), measure.vars = c("w1_gcp", "w2_gcp"), 
                               variable.name = "wave", value.name = "gcp")
    lsim_dt[, wave := as.numeric(str_extract(wave, "[0-9]"))]
    lsim_dt[, time := ifelse(wave == 1, 0, 5)]

    ## run models with and without mortality
    model_base <- lmer(gcp ~ time * any_school + (1 | prim_key), data = lsim_dt)
    params_base <- as.data.table(parameters::model_parameters(model_base))
    model_mort <- lmer(gcp ~ time * any_school + (1 | prim_key), data = lsim_dt[!(died == 1 & wave == 2)])
    params_mort <- as.data.table(parameters::model_parameters(model_mort))

    ## format results
    params_select <- grepl("time:any_school",params_base$Parameter)
    result_dt <- data.table(model = c("Base", "Mortality"),
                            int_value = c(params_base$Coefficient[params_select], params_mort$Coefficient[params_select]), 
                            prop_anyschool = prop_anyschool, cog_anyschool = cog_anyschool, w2gcp_w1gcp = w2gcp_w1gcp, 
                            died_w1gcp = died_w1gcp, died_anyschool = died_anyschool, died_int = died_int, prop_died = num_died/sample_size)
    return(result_dt)
}

## run base simulation 
sim_results <- rbindlist(pbreplicate(iterations, sim_data(), simplify = FALSE))

## extreme effects 
sim_results_extreme <- rbindlist(pbreplicate(iterations, 
                                sim_data(intercept = intercept_extreme, 
                                         gcp_effect = log(0.5),
                                         anyschool_effect = log(0.5),
                                         interaction_effect = log(0.5)),
                                simplify = FALSE))

# PRINT SIMULATION RESULTS -----------------------------------------------------------

## check simulation parameterization
sim_results[model == "Base", lapply(.SD, mean), .SDcols = c("prop_anyschool", "cog_anyschool", "w2gcp_w1gcp", "died_w1gcp", "died_anyschool", "died_int", "prop_died")]
sim_results_extreme[model == "Base", lapply(.SD, mean), .SDcols = c("prop_anyschool", "cog_anyschool", "w2gcp_w1gcp", "died_w1gcp", "died_anyschool", "died_int", "prop_died")]
anyschool_mean; w1cog_params; w2cog_params; survmodel_params

## show and save results
sim_results[, .(mean = mean(int_value), lower = quantile(int_value, 0.025), upper = quantile(int_value, 0.975)), by = "model"]
sim_results_extreme[, .(mean = mean(int_value), lower = quantile(int_value, 0.025), upper = quantile(int_value, 0.975)), by = "model"]

all_results <- rbind(
    sim_results[, .(scenario = "base", mean = mean(int_value), lower = quantile(int_value, 0.025), upper = quantile(int_value, 0.975)), by = "model"], 
    sim_results_extreme[, .(scenario = "extreme", mean = mean(int_value), lower = quantile(int_value, 0.025), upper = quantile(int_value, 0.975)), by = "model"]
    )

write.xlsx(all_results, paste0(plot_dir, "simulation_results_", date, ".xlsx"))
