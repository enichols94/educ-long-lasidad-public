##########################################################################
### Author: Emma Nichols
### Date: 12/11/2024
### Project: LASIDAD Education and longitudinal change
### Purpose: Attrition and death weights
##########################################################################

rm(list = ls())

pacman::p_load(data.table, openxlsx, haven, readr, dplyr, magrittr, stringr,
               labelled, scales, gridExtra, grid, ggplot2, lubridate, 
               survival, patchwork, srvyr, pROC, Hmisc, weights)
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

data <- read_rds(paste0(derived_dir, "processed_data.rds"))
survival_dt <- data$survival
longitudinal_dt <- data$longitudinal

dt <- copy(survival_dt)
dt <- merge(dt, longitudinal_dt[wave == 1, .(prim_key, gcp, memory, executive, language, orientation, visuospatial)], by = "prim_key")

# REFORMAT SOME CONTINUOUS ---------------------------------------------------

dt[, age10 := age/10]

# REDEFINE ATTRITION AS ATTRITION + REFUSAL -------------------------------

dt[, attrited := as.numeric(attrited == 1 | refused == 1)]

# COVS --------------------------------------------------------------------

covs <- c("age10", "female", "statef", "caste", "rural", "educ", "educ_mom", "educ_dad", "wealthq", "consumptionq", "puccahouse", 
          "improvedsanitation", "childhood_finance", "childhood_health",
          "gcp", "memory", "executive", "language", "visuospatial")

# CREATE MISSING INDICATORS -----------------------------------------------

## create missing indicators for any predictors with missingness
for (cov in covs){
  if (nrow(dt[is.na(get(cov))]) > 0){
    if (dt[, class(get(cov))] == "numeric"){
      dt[, paste0(cov, "_ind") := as.numeric(is.na(get(cov)))]
      dt[, paste0(cov, "_reg") := ifelse(is.na(get(cov)), 0, get(cov))]
    } else if (dt[, class(get(cov))] == "factor"){
      dt[, c(cov) := factor(case_when(is.na(get(cov)) ~ "Missing",
                                      !is.na(get(cov)) ~ as.character(get(cov))),
                            levels = c(levels(get(cov)), "Missing"))]
    }
  }
}

## update covariates for model formula
reg_covs <- copy(covs)
for (cov in covs){
  if (paste0(cov, "_reg") %in% names(dt)){
    reg_covs <- reg_covs[!reg_covs == cov]
    reg_covs <- c(reg_covs, paste0(cov, "_reg"), paste0(cov, "_ind"))
  }
}


# MODEL DEATH -----------------------------------------------------------

death_reg_covs <- reg_covs

death_model0 <- glm(as.formula(paste0("died ~", paste(death_reg_covs, collapse = " + "))),
                    data = dt, family = binomial(link = "logit"))
death_auc0 <- auc(roc(response = dt[, died], predictor = predict(death_model0, type = "response")))

## coefficients to test
test_covs <- death_reg_covs[!grepl("_ind$", death_reg_covs)]

## test whether quadratic terms are beneficial
death_quadratic_models <- data.table()
for (cov in test_covs){
  if (dt[, class(get(cov))] == "numeric" & dt[, length(unique(get(cov)))] > 10){
    print(cov)
    death_model_test <- update(death_model0, as.formula(paste0(".~. + I(", cov, "^2)")))
    anova_val <- anova(death_model0, death_model_test, test = "Chisq")

    ## if significant anova - save effect and add to death_reg_covs
    if (anova_val$`Pr(>Chi)`[2] < 0.1){
      death_reg_covs <- c(death_reg_covs, paste0("I(", cov, "^2)"))
      params <- as.data.table(parameters::model_parameters(death_model_test, ci_method = "wald"))
      params[, quadratic_var := cov]
      death_quadratic_models <- rbind(death_quadratic_models, params)
    }
  }
} 

## test whether gender interactions are beneficial
death_gender_models <- data.table()
for (cov in test_covs[!test_covs == "female"]){
  print(cov)
  death_model_test <- update(death_model0, as.formula(paste0(".~. + female:", cov)))
  anova_val <- anova(death_model0, death_model_test, test = "Chisq")

  ## if significant anova - save effect and add to death_reg_covs
  if (anova_val$`Pr(>Chi)`[2] < 0.1){
    death_reg_covs <- c(death_reg_covs, paste0("female:", cov))
    params <- as.data.table(parameters::model_parameters(death_model_test, ci_method = "wald"))
    params[, gender_var := cov]
    death_gender_models <- rbind(death_gender_models, params)
  }
} 

## test whether rural interactions are beneficial
death_rural_models <- data.table()
for (cov in test_covs[!test_covs == "rural"]){
  print(cov)
  death_model_test <- update(death_model0, as.formula(paste0(".~. + rural:", cov)))
  anova_val <- anova(death_model0, death_model_test, test = "Chisq")

  ## if significant anova - save effect and add to death_reg_covs
  if (anova_val$`Pr(>Chi)`[2] < 0.1){
    death_reg_covs <- c(death_reg_covs, paste0("rural:", cov))
    params <- as.data.table(parameters::model_parameters(death_model_test, ci_method = "wald"))
    params[, rural_var := cov]
    death_rural_models <- rbind(death_rural_models, params)
  }
} 

death_model <- glm(as.formula(paste0("died ~", paste(death_reg_covs, collapse = " + "))),
                    data = dt, family = binomial(link = "logit"))
death_auc1 <- auc(roc(response = dt[, died], predictor = predict(death_model, type = "response")))

# MODEL ATTRITION -----------------------------------------------------------

attrition_reg_covs <- reg_covs

attrition_model0 <- glm(as.formula(paste0("attrited ~", paste(attrition_reg_covs, collapse = " + "))),
                    data = dt[died == 0], family = binomial(link = "logit"))
attrition_auc0 <- auc(roc(response = dt[died == 0, attrited], predictor = predict(attrition_model0, type = "response")))

## coefficients to test
test_covs <- attrition_reg_covs[!grepl("_ind$", attrition_reg_covs)]

## test whether quadratic terms are beneficial
attrition_quadratic_models <- data.table()
for (cov in test_covs){
  if (dt[, class(get(cov))] == "numeric" & dt[, length(unique(get(cov)))] > 10){
    print(cov)
    attrition_model_test <- update(attrition_model0, as.formula(paste0(".~. + I(", cov, "^2)")))
    anova_val <- anova(attrition_model0, attrition_model_test, test = "Chisq")

    ## if significant anova - save effect and add to attrition_reg_covs
    if (anova_val$`Pr(>Chi)`[2] < 0.1){
      attrition_reg_covs <- c(attrition_reg_covs, paste0("I(", cov, "^2)"))
      params <- as.data.table(parameters::model_parameters(attrition_model_test, ci_method = "wald"))
      params[, quadratic_var := cov]
      attrition_quadratic_models <- rbind(attrition_quadratic_models, params)
    }
  }
} 

## test whether gender interactions are beneficial
attrition_gender_models <- data.table()
for (cov in test_covs[!test_covs == "female"]){
  print(cov)
  attrition_model_test <- update(attrition_model0, as.formula(paste0(".~. + female:", cov)))
  anova_val <- anova(attrition_model0, attrition_model_test, test = "Chisq")

  ## if significant anova - save effect and add to attrition_reg_covs
  if (anova_val$`Pr(>Chi)`[2] < 0.1){
    attrition_reg_covs <- c(attrition_reg_covs, paste0("female:", cov))
    params <- as.data.table(parameters::model_parameters(attrition_model_test, ci_method = "wald"))
    params[, gender_var := cov]
    attrition_gender_models <- rbind(attrition_gender_models, params)
  }
} 

## test whether rural interactions are beneficial
attrition_rural_models <- data.table()
for (cov in test_covs[!test_covs == "rural"]){
  print(cov)
  attrition_model_test <- update(attrition_model0, as.formula(paste0(".~. + rural:", cov)))
  anova_val <- anova(attrition_model0, attrition_model_test, test = "Chisq")

  ## if significant anova - save effect and add to attrition_reg_covs
  if (anova_val$`Pr(>Chi)`[2] < 0.1){
    attrition_reg_covs <- c(attrition_reg_covs, paste0("rural:", cov))
    params <- as.data.table(parameters::model_parameters(attrition_model_test, ci_method = "wald"))
    params[, rural_var := cov]
    attrition_rural_models <- rbind(attrition_rural_models, params)
  }
} 

attrition_model <- glm(as.formula(paste0("attrited ~", paste(attrition_reg_covs, collapse = " + "))),
                    data = dt[died == 0], family = binomial(link = "logit"))
attrition_auc1 <- auc(roc(response = dt[died == 0, attrited], predictor = predict(attrition_model, type = "response")))

# CALCULATE WEIGHTS -----------------------------------------------------------

dt[, death_num := (died)*dt[, mean(died)] + (1-died)*dt[, 1-mean(died)]]
dt[, death_den := (died)*(predict(death_model, type = "response")) + (1-died)*(1-predict(death_model, type = "response"))]
dt[, death_weight := death_num/death_den]

dt[died == 0, attrition_num := (attrited)*(dt[died == 0, mean(attrited)]) + (1-attrited)*(dt[died == 0, 1-mean(attrited)])]
dt[died == 0, attrition_den := (attrited)*predict(attrition_model, type = "response") + (1-attrited)*(1-predict(attrition_model, type = "response"))]
dt[died == 0, attrition_weight := attrition_num/attrition_den]

dt[died == 0 & attrited == 0, combo_weight := death_weight * attrition_weight]
dt[died == 0 & attrited == 0, full_weight := death_weight * attrition_weight * weight_w1]

dt[, no_weight := 1]

# GRAPHS FOR COVARIATE BALANCE ------------------------------------------------

## use the unweighted pooled variance for the standardized mean difference 
# https://cran.r-project.org/web/packages/cobalt/vignettes/cobalt.html#variance-in-standardized-mean-differences-and-correlations

get_diff_continuous <- function(v, weight = w, data_full = dt_full, data_selected = dt_selected){ ## standardized mean difference for continuous variables
  (weighted.mean(data_selected[!is.na(get(v)), get(v)], w = data_selected[!is.na(get(v)), get(weight)]) - 
     weighted.mean(data_full[!is.na(get(v)), get(v)], w = rep(1, nrow(data_full[!is.na(get(v))]))))/
    (sqrt((wtd.var(data_full[!is.na(get(v)), get(v)], w = rep(1, nrow(data_full[!is.na(get(v))]))) +
             wtd.var(data_selected[!is.na(get(v)), get(v)], w = rep(1, nrow(data_selected[!is.na(get(v))]))))/2))
}

get_diff_binary <- function(v, weight, data_full = dt_full, data_selected = dt_selected){ ## raw mean difference for binary variables
  (weighted.mean(data_selected[!is.na(get(v)), get(v)], w = data_selected[!is.na(get(v)), get(weight)]) - 
     weighted.mean(data_full[!is.na(get(v)), get(v)], w = rep(1, nrow(data_full[!is.na(get(v))]))))
}

get_diff_factor <- function(v, weight, data_full = dt_full, data_selected = dt_selected){ ## raw mean difference for each category of factor variables
  ls <- data_full[, levels(get(v))]
  normal <- wpct(data_full[!is.na(get(v)), get(v)], weight = rep(1, nrow(data_full[!is.na(get(v))])))
  alt <- wpct(data_selected[!is.na(get(v)), get(v)], weight = data_selected[!is.na(get(v)), get(weight)])
  result <- sapply(ls, function(x) unname(alt[names(alt) == x]) - unname(normal[names(normal) == x]))
}

get_meandiff <- function(var, dt_full, dt_selected, w){
  if (!class(dt[, get(var)]) == "factor" & length(unique(dt[!is.na(get(var)), get(var)])) > 2){ ## continuous
    unweighted_diff <- get_diff_continuous(v = var, weight = "no_weight", data_full = dt_full, data_selected = dt_selected)
    weighted_diff <- get_diff_continuous(v = var, weight = w, data_full = dt_full, data_selected = dt_selected)
    result_dt <- data.table(variable = var, diff = c(unweighted_diff, weighted_diff), type = c("Unweighted", "Weighted"))
  } else if (length(unique(dt[!is.na(get(var)), get(var)])) == 2) { ## binary
    unweighted_diff <- get_diff_binary(v = var, weight = "no_weight", data_full = dt_full, data_selected = dt_selected)
    weighted_diff <- get_diff_binary(v = var, weight = w, data_full = dt_full, data_selected = dt_selected)
    result_dt <- data.table(variable = var, diff = c(unweighted_diff, weighted_diff), type = c("Unweighted", "Weighted"))
  } else { ## factors
    ls <- dt[!is.na(get(var)), levels(get(var))]; nlevels <- length(ls)
    unweighted_diff <- get_diff_factor(v = var, weight = "no_weight", data_full = dt_full, data_selected = dt_selected)
    weighted_diff <- get_diff_factor(v = var, weight = w, data_full = dt_full, data_selected = dt_selected)
    
    result_dt <- data.table(variable = var, level = ls, diff = c(unweighted_diff, weighted_diff), 
                            type = c(rep("Unweighted", nlevels), rep("Weighted", nlevels)))
  }
  result_dt[, weight := w]
  return(result_dt)
}

diff_covs <- covs[!covs == "statef"]

death_diffs <- rbindlist(lapply(diff_covs, function(x) get_meandiff(var = x, dt_full = dt, dt_selected = dt[died == 0],
                                                                    w = "death_weight")), fill = T, use.names = T)
attrition_diffs <- rbindlist(lapply(diff_covs, function(x) get_meandiff(var = x, dt_full = dt[died == 0], dt_selected = dt[died == 0 & attrited == 0],
                                                                    w = "attrition_weight")), fill = T, use.names = T)
full_diffs <- rbindlist(lapply(diff_covs, function(x) get_meandiff(var = x, dt_full = dt, dt_selected = dt[died == 0 & attrited == 0],
                                                                    w = "combo_weight")), fill = T, use.names = T)

# CREATE BALANCE GRAPHS -----------------------------------------------------------

balance_graph <- function(diffs, name = ""){
    ind_diffs <- copy(diffs)
    var_conversion <- data.table(variable = c("age10", "female", "caste", "rural", "educ", "educ_mom", "educ_dad", "wealthq", "consumptionq", "puccahouse", 
                                            "improvedsanitation", "childhood_finance", "childhood_health", "gcp", "memory", "executive", "language", "orientation", "visuospatial"),
                                variable_name = c("Age", "Female", "Caste", "Rural", "Education", "Education (maternal)", "Education (paternal)", "Wealth quintile", 
                                                "Consumption quintile", "Formal housing material", "Improved sanitation", "Childhood financial conditions", "Childhood health", 
                                                "General cognitive functioning", "Memory", "Executive functioning", "Language/fluency", "Orientation", "Visuospatial functioning"))
    ind_diffs <- merge(ind_diffs, var_conversion, by = "variable", sort = F)
    ind_diffs[, axis_label := ifelse(is.na(level), variable_name, paste0(variable_name, ": ", level))]
    ind_diffs[, axis_label := factor(axis_label, levels = rev(unique(ind_diffs[, axis_label])))]

    balance_graph <- ggplot(ind_diffs, aes(x = as.factor(axis_label), y = diff, color = as.factor(type))) +
        geom_point(size = 2, alpha = 0.8) + 
        geom_hline(yintercept = 0) +
        geom_hline(yintercept = c(-0.1,0.1), linetype = "dashed") +
        scale_color_manual(name = "", values = c("#103778", "#FF7A48")) +
        coord_flip() +
        labs(y = "Difference", x = "", title = name) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
    return(balance_graph)
}

death_balance <- balance_graph(death_diffs, name = "Death")
attrition_balance <- balance_graph(attrition_diffs, name = "Attrition")
combo_balance <- balance_graph(full_diffs, name = "Death + Attrition")

pdf(paste0(plot_dir, "balance_plots_", date, ".pdf"), height = 8, width = 8)
death_balance; attrition_balance; combo_balance
dev.off()

# ADD WEIGHT TO DATA -----------------------------------------------------------

weight_dt <- copy(dt[, .(prim_key, death_weight, attrition_weight, combo_weight, full_weight)])

survival_dt <- merge(survival_dt, weight_dt, by = "prim_key")
longitudinal_dt <- merge(longitudinal_dt, weight_dt, by = "prim_key")

write_rds(list(survival = survival_dt, longitudinal = longitudinal_dt), 
          paste0(derived_dir, "processed_data_weights.rds"))
