##########################################################################
### Author: Emma Nichols
### Date: 12/11/2024
### Project: LASIDAD Educationa and longitudinal change
### Purpose: Mixed effects models
##########################################################################

rm(list = ls())

pacman::p_load(data.table, openxlsx, haven, readr, dplyr, magrittr, stringr,
               labelled, scales, gridExtra, grid, ggplot2, lubridate, INLAjoint, INLA,
               survival, patchwork, survey, srvyr, lme4, lqmm)
date <- gsub("-", "_", Sys.Date())
set.seed(6541)

## add gender to stratifications

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

model_covs <- c("age", "gender", "caste", "rural", "educ_dad", "childhood_finance", "childhood_health")

# FORMAT DATA -----------------------------------------------------------

## note: maybe should redo weights with this exclusion or create a weight to deal with this

## remove missing data
dt <- copy(longitudinal_dt); surv_dt <- copy(survival_dt)
for (var in c("gcp", "educ", model_covs)){
    message(paste0("Missing individuals in ", var, ": ", dt[is.na(get(var)), length(unique(prim_key))], " (", 
                   dt[is.na(get(var)) & wave == 2, length(unique(prim_key))], " longitudinal)"))
    dt <- dt[!is.na(get(var))]
}
surv_dt <- surv_dt[prim_key %in% dt[, unique(prim_key)]] ## only keep those in longitudinal data

## collapse childhood health
dt[childhood_health == "Very poor", childhood_health := "Poor"]
dt[, childhood_health := droplevels(childhood_health)]
surv_dt[childhood_health == "Very poor", childhood_health := "Poor"]
surv_dt[, childhood_health := droplevels(childhood_health)]

# MODELS -----------------------------------------------------------

## mixed effects model (base)
m0_form <- paste0("gcp ~ educ + educ*time + ", paste0(model_covs, collapse = "*time + "), "*time + (1 | prim_key)")
m0 <- lmer(as.formula(m0_form), data = dt)

## weighted gee models
m1_form <- paste0("gcp ~ educ + educ*time + ", paste0(model_covs, collapse = "*time + "), "*time")
m1 <- geepack::geeglm(as.formula(m1_form), data = dt, id = dt[, prim_key], family = gaussian, weights = dt[, weight_w1*death_weight])
m2 <- geepack::geeglm(as.formula(m1_form), data = dt, id = dt[, prim_key], family = gaussian, weights = dt[, weight_w1*attrition_weight])
m3 <- geepack::geeglm(as.formula(m1_form), data = dt, id = dt[, prim_key], family = gaussian, weights = dt[, full_weight])

## joint model
message("Joint model further excludes ", length(surv_dt[attrited == 1]), " individuals")
m4_form_l <- gsub("prim_key", "id", m0_form)
m4_form_s <- paste0("inla.surv(time_elapsed, died) ~ educ + ", paste0(model_covs, collapse = " + "))
joint_dt <- copy(dt[order(prim_key, wave) & !prim_key %in% surv_dt[attrited == 1, unique(prim_key)]]) ## for joint models get data, sort, and create new ID
joint_surv_dt <- copy(surv_dt[order(prim_key) & !attrited == 1])
joint_surv_dt[, id := 1:.N]
joint_dt <- merge(joint_dt, joint_surv_dt[, .(prim_key, id)], by = "prim_key")    
m4 <- joint(formSurv = as.formula(m4_form_s), formLong = as.formula(m4_form_l), family = "normal", 
            dataLong = joint_dt, dataSurv = joint_surv_dt, 
            id = "id", timeVar = "time_elapsed", assoc = "CV", 
            basRisk  = "rw2", NbasRisk=10, control=list(int.strategy="eb"))

## stratified models - by age group, rural/urban, above/below median W1 cog, above/below median follow-up
get_model <- function(model_data){
    return(lmer(as.formula(m0_form), data = model_data))
}

m5 <- get_model(dt[age_group == "60-69"]); m6 <- get_model(dt[age_group == "70-79"]); m7 <- get_model(dt[age_group == "80+"])
m8 <- get_model(dt[rural == 0]); m9 <- get_model(dt[rural == 1])
m10 <- get_model(dt[female == 0]); m11 <- get_model(dt[female == 1])
m12 <- get_model(dt[above_w1medcog == 0]); m13 <- get_model(dt[above_w1medcog == 1])
m14 <- get_model(dt[above_w2medfu == 0]); m15 <- get_model(dt[above_w2medfu == 1])
m16 <- get_model(dt[cdr == "Normal"]); m17 <- get_model(dt[cdr == "MCI/Questionnable"]); m18 <- get_model(dt[cdr == "Dementia"])

## quantile regression
m19 <- rq(as.formula(m1_form), data = dt, tau = c(0.1,0.25,0.5,0.75, 0.9))
m19_params <- summary(m19, se = "boot", cluster = dt[, prim_key]) ## https://ideas.repec.org/a/taf/jnlasa/v112y2017i517p446-456.html

## individual test items
get_model_item <- function(outcome_item){
    return(lmer(as.formula(gsub("gcp", outcome_item, m0_form)), data = dt))
}

m20 <- get_model_item("memory"); m21 <- get_model_item("executive"); m22 <- get_model_item("language"); m23 <- get_model_item("visuospatial")
m24 <- get_model_item("wr_delayed"); m25 <- get_model_item("lm_delayed"); m26 <- get_model_item("animals"); m27 <- get_model_item("con_praxis"); m28 <- get_model_item("ravens")

# FORMAT MODEL RESULTS -----------------------------------------------------------

get_coefficients <- function(model_num){
    model <- get(paste0("m", model_num))

    if (model_num == 4){ ## alternative for joint model (m4)
        params <- summary(model)$FixedEff[[1]]
        param_select <- grepl("^educ.*\\:time_L1$", rownames(params))
        result_dt <- data.table(model_num = model_num, 
                                parameter = rownames(params)[param_select], 
                                estimate = params[param_select, "mean"], 
                                lower = params[param_select, "0.025quant"], 
                                upper = params[param_select, "0.975quant"])
    } else if (model_num == 19){ ## alternative for quantile regression (m19)
        params <- m19_params

        get_coef <- function(tau_num){
            selected <- params[[tau_num]]
            param_select <- grepl("^educ.*\\:time$", row.names(selected$coefficients))
            select_dt <- data.table(model_num = paste0(model_num, " - ", selected$tau), 
                                    parameter = row.names(selected$coefficients)[param_select], 
                                    estimate = data.frame(selected$coefficients)$Value[param_select], 
                                    lower = data.frame(selected$coefficients)$Value[param_select]-1.96*data.frame(selected$coefficients)$`Std..Error`[param_select], ## assume CLT, CI symmetry
                                    upper = data.frame(selected$coefficients)$Value[param_select]+1.96*data.frame(selected$coefficients)$`Std..Error`[param_select])
            return(select_dt)
        }

        result_dt <- rbindlist(lapply(1:length(params), get_coef))


    } else { ## all other models
        params <- parameters::parameters(model)
        param_select <- grepl("^educ.*\\:time$", params$Parameter)
        result_dt <- data.table(model_num = model_num, 
                                parameter = params$Parameter[param_select], 
                                estimate = params$Coefficient[param_select], 
                                lower = params$CI_low[param_select], 
                                upper = params$CI_high[param_select])
    }
    
    return(result_dt)
}

model_results <- rbindlist(lapply(0:28, get_coefficients))

# FORMAT GRAPH DATA -----------------------------------------------------------

## make parameter names consistent
model_results[, parameter := gsub("_L1$", "", parameter)][, parameter := gsub("\\.", " ", parameter)][parameter == "educMiddlesecondary school:time", parameter := "educMiddle-secondary school:time"] 

## get education categories
model_results[, educ_cat := gsub("educ", "", parameter)][, educ_cat := gsub(":time", "", educ_cat)][, educ_cat := factor(educ_cat, levels = c("Less than primary", "Primary school", "Middle-secondary school", "Higher secondary school and up"))]

## merge on labels and categories
label_dt <- data.table(model_num = model_results[, unique(model_num)], 
                       labels = c("Mixed effects model (base)", "Weighted GEE model (survey + death)", "Weighted GEE model (survey + attrition)", 
                                  "Weighted GEE model (survey + death + attrition)", "Joint model", "Age (60-69)", "Age (70-79)", 
                                  "Age (80+)", "Urbanicity (rural)", "Urbanicity (urban)", "Gender (men)", "Gender (women)",
                                  "W1 cognition (below median W1 cog)", "W1 cognition (above median W1 cog)", 
                                  "W2 follow-up (below median follow-up)", "W2 follow-up (above median follow-up)", 
                                  "W1 CDR (normal)", "W1 CDR (MCI/Questionable)", "W1 CDR (dementia)",
                                  "10th percentile", "25th percentile", "50th percentile", "75th percentile", "90th percentile",
                                  "Memory", "Executive functioning", "Language", "Visuospatial functioning", 
                                  "Delayed word recall", "Delayed story recall (logical memory)", "Animal naming", "Constructional praxis", "Raven's progressive matrices"), 
                       category = c(rep("Base models",5), rep("Stratified", 14), rep("Quantiles", 5), rep("Other cognitive outcomes", 9)))
label_dt[, labels := factor(labels, levels = rev(label_dt[, labels]))]                       
model_results <- merge(model_results, label_dt, by = "model_num", sort = FALSE)

# MAKE GRAPH -----------------------------------------------------------

get_graphs <- function(c){
    plot <- ggplot(model_results[category == c], aes(x = labels, y = estimate, ymin = lower, ymax = upper, color = educ_cat)) +
        geom_point() +
        geom_errorbar(width = 0.1) + 
        geom_hline(yintercept = 0, linetype = "dashed") +
        facet_wrap(~educ_cat, nrow = 1) + 
        coord_flip() + 
        labs(x = "", y = "") +
        scale_y_continuous(limits = c(-0.102,0), breaks = seq(-.1,0,0.025), oob = oob_keep) +
        scale_color_manual(name = "", values = c("#6551CC", "#CC5151", "#2C85B2", "#85B22C")) + 
        theme_bw() +
        theme(legend.position = "none", plot.margin = unit(c(0,5.5,0,5.5), "pt"))
    if (c == "Other cognitive outcomes"){
        plot <- plot + labs(y = "Difference in rate of annual cognitive decline (Ref: No schooling)") + theme(legend.position = "bottom")
    }
    if (!c == "Base models"){
        plot <- plot + theme(strip.text = element_blank())
    }
    if (!c == "Other cognitive outcomes"){
        plot <- plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
    }
    return(plot)
}

plots <- lapply(model_results[, unique(category)], get_graphs)  

get_labels <- function(category){
    label_plot <- ggplot() + 
        theme_void() +
        geom_text(aes(0,0, label = str_wrap(category, width = 12)))
    return(label_plot)
}

labels <- lapply(model_results[, unique(category)], get_labels)

## combine plots
full_plot <- plots[[1]] + labels[[1]] + plots[[2]] + labels[[2]] + plots[[3]] + labels[[3]] + plots[[4]] + labels[[4]] +
    plot_layout(nrow = 4, widths = c(1, 0.1), heights = c(1,2.3,0.9,1.5))

ggsave(paste0(plot_dir, "regression_comparison_", date, ".pdf"), plot= full_plot, height = 8, width = 15)
