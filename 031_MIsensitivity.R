# Create multiple imputations for the longitudinal-model sensitivity analysis.
#
# Input: data/derived/processed_data_weights.rds from 011_weights.R.
# Output: data/derived/mi_data.rds, consumed by 030_longmodels.R.

source(file.path("R", "project_paths.R"))
load_required_packages(c(
  "data.table", "openxlsx", "haven", "readr", "dplyr", "magrittr",
  "stringr", "labelled", "scales", "gridExtra", "grid", "ggplot2",
  "lubridate", "INLAjoint", "INLA", "survival", "patchwork", "survey",
  "srvyr", "lme4", "lqmm", "broom.mixed", "mice"
))
set.seed(6541)

# PATHS -------------------------------------------------------------------

paths <- project_paths()
derived_dir <- paths$derived
ensure_output_directory(derived_dir)

# READ DATA ------------------------------------------------------------------

data <- read_rds(required_input_file(file.path(derived_dir, "processed_data_weights.rds")))
survival_dt <- data$survival
longitudinal_dt <- data$longitudinal

model_covs <- c("age", "female", "caste", "rural", "educ_dad", "childhood_finance", "childhood_health")

# FORMAT DATA -----------------------------------------------------------

## number lost to mortality and loss to follow-up  
survival_dt[, sum(died)]; survival_dt[, sum(as.numeric(refused == 1 | attrited == 1))]

## rename data
dt <- copy(longitudinal_dt)

## collapse childhood health
dt[childhood_health == "Very poor", childhood_health := "Poor"]
dt[, childhood_health := droplevels(childhood_health)]

# SET UP MULTIPLE IMPUTATION ---------------------------------------------

mi_dt <- dcast.data.table(dt, as.formula(paste0("prim_key + educ + ", paste(model_covs[!model_covs == "age"], collapse = " + "), " ~ wave")), 
                          value.var = c("gcp", "age", "time"))

## make sure all binary vars are factors
for (var in names(mi_dt)){
    if (mi_dt[!is.na(get(var)), length(unique(get(var)))] == 2){
        mi_dt[, c(var) := as.factor(get(var))]
    }
}

## create predictor matrix
get_predmat <- function(dt){

    ismissing <- sapply(names(dt), function(x) as.numeric(nrow(dt[is.na(get(x))]) > 0)) ## what variables have missingness? 
    impute_vars <- names(ismissing)[ismissing == 1] ## variables to impute
    impute_vars <- impute_vars[!impute_vars %in% c("gcp_2", "time_2")] ## don't impute these - missingness from missing follow-up

    predictor_matrix <- matrix(0, ncol = length(names(dt)), nrow = length(names(dt))) ## create empty matrix
    predictor_matrix[match(impute_vars, names(dt)),2:(ncol(dt))] <- 1 ## use all variables except the ID as candidate predictors
    predictor_matrix[, match(c("gcp_2", "time_1", "time_2"), names(dt))] <- 0 ## don't use gcp2 or time as a predictor
    diag(predictor_matrix) <- 0 ## prevent variables from predicting themselves

    return(predictor_matrix)
}
pred_matrix <- get_predmat(mi_dt)

## create string of methods - unimputed variables must be empty ""
get_methods <- function(dt){

    ismissing <- sapply(names(dt), function(x) as.numeric(nrow(dt[is.na(get(x))]) > 0)) ## what variables have missingness? 
    impute_vars <- names(ismissing)[ismissing == 1] ## variables to impute
    impute_vars <- impute_vars[!impute_vars %in% c("gcp_2", "time_2")] ## don't impute these
    methods <- c()
    for (var in names(dt)){
        if (var %in% impute_vars){
            if (dt[, class(get(var))] == "numeric"){
                methods <- c(methods, "pmm")
            } else if (dt[, class(get(var))] == "factor" & dt[, length(levels(get(var)))] == 2){
                methods <- c(methods, "logreg")
            } else if (dt[, class(get(var))] == "factor" & dt[, length(levels(get(var)))] > 2){
                methods <- c(methods, "polyreg")
            }
        } else {
            methods <- c(methods, "")
        }
    }
    return(methods)
}
imp_methods <- get_methods(mi_dt)

## run imputations
full_imputations <- mice::mice(mi_dt, m = 10, maxit = 10, predictorMatrix = pred_matrix, method = imp_methods, seed = 8729)
write_rds(full_imputations, file.path(derived_dir, "mi_data.rds"))

plot(full_imputations, layout = c(8,4))

