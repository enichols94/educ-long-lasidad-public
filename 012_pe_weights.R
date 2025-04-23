##########################################################################
### Author: Emma Nichols
### Date: 12/11/2024
### Project: LASIDAD Education and longitudinal change
### Purpose: Practice effects weights
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
plot_dir <- paste0(dir, "paper/appendix_figs/")

# READ DATA ------------------------------------------------------------------

data <- read_rds(paste0(derived_dir, "processed_data.rds"))
pe_dt <- data$practice

# REFORMAT IF NEEDED ---------------------------------------------------------

pe_dt[, age10 := age/10]
pe_dt[, practice := as.numeric(!refresher == 1)]

# COVS -----------------------------------------------------------------------

covs <- c("age10", "female", "statef", "caste", "rural", "educ", "educ_mom", "educ_dad", "wealthq", "consumptionq", "puccahouse", 
          "improvedsanitation", "childhood_finance", "childhood_health",
          "cidisymptoms", "eversmoke", "nowsmoke", "proxy", "adls", "iadls", "gripstrength", "gcp_lasi", 
          "jorm", "blessed2", "blessed1")  

# CREATE MISSING INDICATORS -----------------------------------------------

## create missing indicators for any predictors with missingness
for (cov in covs){
  if (nrow(pe_dt[is.na(get(cov))]) > 0){
    if (pe_dt[, class(get(cov))] == "numeric"){
      pe_dt[, paste0(cov, "_ind") := as.numeric(is.na(get(cov)))]
      pe_dt[, paste0(cov, "_reg") := ifelse(is.na(get(cov)), 0, get(cov))]
    } else if (pe_dt[, class(get(cov))] == "factor"){
      pe_dt[, c(cov) := factor(case_when(is.na(get(cov)) ~ "Missing",
                                      !is.na(get(cov)) ~ as.character(get(cov))),
                            levels = c(levels(get(cov)), "Missing"))]
    }
  }
}

## update covariates for model formula
reg_covs <- copy(covs)
for (cov in covs){
  if (paste0(cov, "_reg") %in% names(pe_dt)){
    reg_covs <- reg_covs[!reg_covs == cov]
    reg_covs <- c(reg_covs, paste0(cov, "_reg"), paste0(cov, "_ind"))
  }
}

## coefficients to test
test_covs <- reg_covs[!grepl("_ind$", reg_covs)]

# MODEL EXPOSURE TO PRACTICE (WHEHTER IN REFRESHER/RETURNER SAMPLE) --------

## reg covs for this round of modeling
reg_covs1 <- copy(reg_covs)

model0 <- glm(as.formula(paste0("practice ~", paste(reg_covs1, collapse = " + "))),
                    data = pe_dt, family = binomial(link = "logit"))
auc0 <- auc(roc(response = pe_dt[, practice], predictor = predict(model0, type = "response")))

## test whether quadratic terms are beneficial
quadratic_models <- data.table()
for (cov in test_covs){
  if (pe_dt[, class(get(cov))] == "numeric" & pe_dt[, length(unique(get(cov)))] > 10){
    print(cov)
    model_test <- update(model0, as.formula(paste0(".~. + I(", cov, "^2)")))
    anova_val <- anova(model0, model_test, test = "Chisq")

    ## if significant anova - save effect and add to reg_covs1
    if (anova_val$`Pr(>Chi)`[2] < 0.1){
      reg_covs1 <- c(reg_covs1, paste0("I(", cov, "^2)"))
      params <- as.data.table(parameters::model_parameters(model_test, ci_method = "wald"))
      params[, quadratic_var := cov]
      quadratic_models <- rbind(quadratic_models, params)
    }
  }
} 

## test whether gender interactions are beneficial
gender_models <- data.table()
for (cov in test_covs[!test_covs == "female"]){
  print(cov)
  model_test <- update(model0, as.formula(paste0(".~. + female:", cov)))
  anova_val <- anova(model0, model_test, test = "Chisq")

  ## if significant anova - save effect and add to reg_covs1
  if (anova_val$`Pr(>Chi)`[2] < 0.1){
    reg_covs1 <- c(reg_covs1, paste0("female:", cov))
    params <- as.data.table(parameters::model_parameters(model_test, ci_method = "wald"))
    params[, gender_var := cov]
    gender_models <- rbind(gender_models, params)
  }
} 

## test whether rural interactions are beneficial
rural_models <- data.table()
for (cov in test_covs[!test_covs == "rural"]){
  print(cov)
  model_test <- update(model0, as.formula(paste0(".~. + rural:", cov)))
  anova_val <- anova(model0, model_test, test = "Chisq")

  ## if significant anova - save effect and add to reg_covs1
  if (anova_val$`Pr(>Chi)`[2] < 0.1){
    reg_covs1 <- c(reg_covs1, paste0("rural:", cov))
    params <- as.data.table(parameters::model_parameters(model_test, ci_method = "wald"))
    params[, rural_var := cov]
    rural_models <- rbind(rural_models, params)
  }
} 

model1 <- glm(as.formula(paste0("practice ~", paste(reg_covs1, collapse = " + "))),
                    data = pe_dt, family = binomial(link = "logit"))
auc1 <- auc(roc(response = pe_dt[, practice], predictor = predict(model1, type = "response")))

# CALCULATE WEIGHTS -----------------------------------------------------------

## propensity score
pe_dt[, ps := predict(model1, type = "response")]

pe_dt[, pe_num := (practice)*pe_dt[, mean(practice)] + (1-practice)*pe_dt[, 1-mean(practice)]]
pe_dt[, pe_den := (practice)*(predict(model1, type = "response")) + (1-practice)*(1-predict(model1, type = "response"))]
pe_dt[, pe_weight := pe_num/pe_den]

pe_dt[, no_weight := 1]

# GRAPHS FOR COVARIATE BALANCE ------------------------------------------------

## use the unweighted pooled variance for the standardized mean difference 
# https://cran.r-project.org/web/packages/cobalt/vignettes/cobalt.html#variance-in-standardized-mean-differences-and-correlations

get_diff_continuous <- function(v, weight = w, data = dt){ ## standardized mean difference for continuous variables
  (weighted.mean(data[!is.na(get(v)) & practice == 1, get(v)], w = data[!is.na(get(v)) & practice == 1, get(weight)]) - 
     weighted.mean(data[!is.na(get(v)) & practice == 0, get(v)], w = data[!is.na(get(v)) & practice == 0, get(weight)]))/
    (sqrt((wtd.var(data[!is.na(get(v)) & practice == 1, get(v)], w = rep(1, nrow(data[!is.na(get(v)) & practice == 1]))) +
             wtd.var(data[!is.na(get(v)) & practice == 0, get(v)], w = rep(1, nrow(data[!is.na(get(v)) & practice == 0]))))/2))
}

get_diff_binary <- function(v, weight, data = dt){ ## raw mean difference for binary variables
  (weighted.mean(data[!is.na(get(v)) & practice == 1, get(v)], w = data[!is.na(get(v)) & practice == 1, get(weight)]) - 
     weighted.mean(data[!is.na(get(v)) & practice == 0, get(v)], w = data[!is.na(get(v)) & practice == 0, get(weight)]))
}

get_diff_factor <- function(v, weight, data = dt){ ## raw mean difference for each category of factor variables
  ls <- data[, levels(get(v))]
  normal <- wpct(data[!is.na(get(v)) & practice == 0, get(v)], weight = data[!is.na(get(v)) & practice == 0, get(weight)])
  alt <- wpct(data[!is.na(get(v)) & practice == 1, get(v)], weight = data[!is.na(get(v)) & practice == 1, get(weight)])
  result <- sapply(ls, function(x) unname(alt[names(alt) == x]) - unname(normal[names(normal) == x]))
}

get_meandiff <- function(var, dt, w){
  if (!class(dt[, get(var)]) == "factor" & length(unique(dt[!is.na(get(var)), get(var)])) > 2){ ## continuous
    unweighted_diff <- get_diff_continuous(v = var, weight = "no_weight", data = dt)
    weighted_diff <- get_diff_continuous(v = var, weight = w, data = dt)
    result_dt <- data.table(variable = var, diff = c(unweighted_diff, weighted_diff), type = c("Unweighted", "Weighted"))
  } else if (length(unique(dt[!is.na(get(var)), get(var)])) == 2) { ## binary
    unweighted_diff <- get_diff_binary(v = var, weight = "no_weight", data = dt)
    weighted_diff <- get_diff_binary(v = var, weight = w, data = dt)
    result_dt <- data.table(variable = var, diff = c(unweighted_diff, weighted_diff), type = c("Unweighted", "Weighted"))
  } else { ## factors
    ls <- dt[!is.na(get(var)), levels(get(var))]; nlevels <- length(ls)
    unweighted_diff <- get_diff_factor(v = var, weight = "no_weight", data = dt)
    weighted_diff <- get_diff_factor(v = var, weight = w, data = dt)
    
    result_dt <- data.table(variable = var, level = ls, diff = c(unweighted_diff, weighted_diff), 
                            type = c(rep("Unweighted", nlevels), rep("Weighted", nlevels)))
  }
  result_dt[, weight := w][, diff := as.numeric(diff)]
  return(result_dt)
}

diff_covs <- covs[!covs == "statef"]

pe_diffs <- rbindlist(lapply(diff_covs, function(x) get_meandiff(var = x, dt = pe_dt, w = "pe_weight")), fill = T, use.names = T)

# CREATE BALANCE GRAPHS -----------------------------------------------------------

balance_graph <- function(diffs, name = ""){
    ind_diffs <- copy(diffs)
    var_conversion <- data.table(variable = c("age10", "female", "caste", "rural", "educ", "educ_mom", "educ_dad", "wealthq", "consumptionq", "puccahouse", 
                                              "improvedsanitation", "childhood_finance", "childhood_health",
                                              "cidisymptoms", "eversmoke", "nowsmoke", "proxy", "adls", "iadls", "gripstrength", "gcp_lasi", 
                                              "jorm", "blessed2", "blessed1"),
                                variable_name = c("Age", "Female", "Caste", "Rural", "Education", "Education (maternal)", "Education (paternal)", "Wealth quintile", "Consumption quintile", "Formal housing material",
                                                  "Improved sanitation", "Childhood financial conditions", "Childhood health", 
                                                  "Depressive symptoms", "Ever smoking", "Current smoking", "Proxy respondent in LASI", "ADLs", "IADLs", "Grip strength", "Cognition", 
                                                  "Jorm (LASI-DAD)", "Blessed Part 2 (LASI-DAD)", "Blessed Part 1 (LASI-DAD)"))
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

pe_balance <- balance_graph(pe_diffs, name = "Practice effects")

# LOOK AT PROPENSITY SCORES ---------------------------------------------------

density_gg <- ggplot(pe_dt, aes(x = ps, fill = as.factor(practice))) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(name = "", values = c("#04BF8A", "#026873"), labels = c("0" = "Refresher", "1" = "Returner")) +
    labs(x = "Propensity score", y = "Density") +
    theme_bw() +
    theme(legend.position = "bottom")

# TRY WITH RESTRICTION --------------------------------------------------------

## implement Strummer trimming and then re-estimate propensity scores afterwards
## https://pmc.ncbi.nlm.nih.gov/articles/PMC8327194/

## 5th percentile in the treated
lower_cut <- pe_dt[practice == 0, quantile(ps, probs = 0.05, na.rm = T)]

## 95th percentile in the untreated
upper_cut <- pe_dt[practice == 1, quantile(ps, probs = 0.95, na.rm = T)]

## restrict and created trimmed dataset
pe_restricted_dt <- copy(pe_dt[ps > lower_cut & ps < upper_cut])

## number trimmed and sample size 
nrow(pe_restricted_dt); nrow(pe_dt)-nrow(pe_restricted_dt)

# CALCULATE NEW WEIGHTS ---------------------------------------------------------

## reg covs for this round of modeling
reg_covs2 <- copy(reg_covs)

model2 <- glm(as.formula(paste0("practice ~", paste(reg_covs2, collapse = " + "))),
                    data = pe_restricted_dt, family = binomial(link = "logit"))
auc2 <- auc(roc(response = pe_restricted_dt[, practice], predictor = predict(model2, type = "response")))

## test whether quadratic terms are beneficial
quadratic_models <- data.table()
for (cov in test_covs){
  if (pe_restricted_dt[, class(get(cov))] == "numeric" & pe_restricted_dt[, length(unique(get(cov)))] > 10){
    print(cov)
    model_test <- update(model2, as.formula(paste0(".~. + I(", cov, "^2)")))
    anova_val <- anova(model2, model_test, test = "Chisq")

    ## if significant anova - save effect and add to reg_covs2
    if (anova_val$`Pr(>Chi)`[2] < 0.1){
      reg_covs2 <- c(reg_covs2, paste0("I(", cov, "^2)"))
      params <- as.data.table(parameters::model_parameters(model_test, ci_method = "wald"))
      params[, quadratic_var := cov]
      quadratic_models <- rbind(quadratic_models, params)
    }
  }
} 

## test whether gender interactions are beneficial
gender_models <- data.table()
for (cov in test_covs[!test_covs == "female"]){
  print(cov)
  model_test <- update(model2, as.formula(paste0(".~. + female:", cov)))
  anova_val <- anova(model2, model_test, test = "Chisq")

  ## if significant anova - save effect and add to reg_covs2
  if (anova_val$`Pr(>Chi)`[2] < 0.1){
    reg_covs2 <- c(reg_covs2, paste0("female:", cov))
    params <- as.data.table(parameters::model_parameters(model_test, ci_method = "wald"))
    params[, gender_var := cov]
    gender_models <- rbind(gender_models, params)
  }
} 

## test whether rural interactions are beneficial
rural_models <- data.table()
for (cov in test_covs[!test_covs == "rural"]){
  print(cov)
  model_test <- update(model2, as.formula(paste0(".~. + rural:", cov)))
  anova_val <- anova(model2, model_test, test = "Chisq")

  ## if significant anova - save effect and add to reg_covs2
  if (anova_val$`Pr(>Chi)`[2] < 0.1){
    reg_covs2 <- c(reg_covs2, paste0("rural:", cov))
    params <- as.data.table(parameters::model_parameters(model_test, ci_method = "wald"))
    params[, rural_var := cov]
    rural_models <- rbind(rural_models, params)
  }
} 

model3 <- glm(as.formula(paste0("practice ~", paste(reg_covs2, collapse = " + "))),
                    data = pe_restricted_dt, family = binomial(link = "logit"))
auc3 <- auc(roc(response = pe_restricted_dt[, practice], predictor = predict(model3, type = "response")))

# EVALUATE/CHECK -------------------------------------------------------------

## calculate propensity score/weights
pe_restricted_dt[, ps := predict(model3, type = "response")]
pe_restricted_dt[, pe_num := (practice)*pe_restricted_dt[, mean(practice)] + (1-practice)*pe_restricted_dt[, 1-mean(practice)]]
pe_restricted_dt[, pe_den := (practice)*(predict(model3, type = "response")) + (1-practice)*(1-predict(model3, type = "response"))]
pe_restricted_dt[, pe_weight := pe_num/pe_den]
pe_restricted_dt[, no_weight := 1]

pe_restricted_diffs <- rbindlist(lapply(diff_covs, function(x) get_meandiff(var = x, dt = pe_restricted_dt, w = "pe_weight")), fill = T, use.names = T)
pe_restricted_balance <- balance_graph(pe_restricted_diffs, name = "Practice effects (restricted sample)")

## density graph after restriction
density_restricted_gg <- ggplot(pe_restricted_dt, aes(x = ps, fill = as.factor(practice))) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(name = "", values = c("#04BF8A", "#026873"), labels = c("0" = "Refresher", "1" = "Returner")) +
    labs(x = "Propensity score", y = "Density") +
    theme_bw() +
    theme(legend.position = "bottom")

pdf(paste0(plot_dir, "practiceeffects_balance_plots_", date, ".pdf"), height = 8, width = 14)
pe_balance; pe_restricted_balance
dev.off()

pdf(paste0(plot_dir, "practiceeffects_density_plots_", date, ".pdf"), height = 4, width = 4)
density_gg; density_restricted_gg
dev.off()

# SAVE DATA -----------------------------------------------------------

write_rds(pe_restricted_dt, 
          paste0(derived_dir, "processed_practiceeffects_data.rds"))
