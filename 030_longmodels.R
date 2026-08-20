# Fit the primary longitudinal models and specified sensitivity analyses.
#
# Inputs: processed_data_weights.rds from 011_weights.R and mi_data.rds from
# 031_MIsensitivity.R. Outputs are date-stamped model figures and tables.

source(file.path("R", "project_paths.R"))
load_required_packages(c(
  "data.table", "openxlsx", "haven", "readr", "dplyr", "magrittr",
  "stringr", "geepack", "labelled", "scales", "gridExtra", "grid",
  "ggplot2", "lubridate", "INLAjoint", "INLA", "survival", "patchwork",
  "survey", "srvyr", "lme4", "lqmm", "mice", "miceadds", "broom.mixed",
  "parameters", "quantreg", "forcats"
))
date <- gsub("-", "_", Sys.Date())
set.seed(6541)

# PATHS -------------------------------------------------------------------

paths <- project_paths()
derived_dir <- paths$derived
plot_dir <- paths$model_figures
appendix_dir <- paths$appendix
ensure_output_directory(plot_dir)
ensure_output_directory(appendix_dir)

# READ DATA ------------------------------------------------------------------

data <- read_rds(required_input_file(file.path(derived_dir, "processed_data_weights.rds")))
survival_dt <- data$survival
longitudinal_dt <- data$longitudinal

model_covs <- c("age", "female", "caste", "rural", "educ_dad", "childhood_finance", "childhood_health")

# FORMAT DATA -----------------------------------------------------------

## number lost to mortality and loss to follow-up  
survival_dt[, sum(died)]; survival_dt[, sum(as.numeric(refused == 1 | attrited == 1))]

## remove missing data
dt <- copy(longitudinal_dt); surv_dt <- copy(survival_dt)
for (var in c("gcp", "educ", model_covs)){
    message(paste0("Missing individuals in ", var, ": ", dt[is.na(get(var)), length(unique(prim_key))], " (", 
                   dt[is.na(get(var)) & wave == 2, length(unique(prim_key))], " longitudinal)"))
    dt <- dt[!is.na(get(var))]
    message(paste0("Now excluding N=", nrow(surv_dt[!prim_key %in% dt[, unique(prim_key)]]), " from data completely"))
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

## cognitive level results 
parameters::parameters(m0)

## weighted gee models
m1_form <- paste0("gcp ~ educ + educ*time + ", paste0(model_covs, collapse = "*time + "), "*time")
m1 <- geepack::geeglm(as.formula(m1_form), data = dt, id = dt[, prim_key], family = gaussian, corstr = "exchangeable", weights = dt[, death_weight]) ## only include death weight
m2 <- geepack::geeglm(as.formula(m1_form), data = dt, id = dt[, prim_key], family = gaussian, corstr = "exchangeable", weights = dt[, weight_w1*attrition_weight])
m3 <- geepack::geeglm(as.formula(m1_form), data = dt, id = dt[, prim_key], family = gaussian, corstr = "exchangeable", weights = dt[, full_weight])

## joint model
message("Joint model further excludes ", length(surv_dt[attrited == 1]), " individuals")
m4_form_l <- gsub("prim_key", "id", m0_form)
m4_form_s <- paste0("inla.surv(time_elapsed, died) ~ educ + ", paste0(model_covs, collapse = " + "))
joint_dt <- copy(dt[order(prim_key, wave) & !prim_key %in% surv_dt[attrited == 1, unique(prim_key)]]) ## for joint models get data, sort, and create new ID
joint_surv_dt <- copy(surv_dt[order(prim_key) & !attrited == 1])
joint_surv_dt[, id := 1:.N]
joint_dt <- merge(joint_dt, joint_surv_dt[, .(prim_key, id)], by = "prim_key")    
m4 <- joint(formSurv = as.formula(m4_form_s), formLong = as.formula(m4_form_l), family = "normal", 
            dataLong = as.data.frame(joint_dt), dataSurv = as.data.frame(joint_surv_dt), 
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
m19 <- quantreg::rq(as.formula(m1_form), data = dt, tau = c(0.1,0.25,0.5,0.75, 0.9))
m19_params <- summary(m19, se = "boot", cluster = dt[, prim_key]) ## https://ideas.repec.org/a/taf/jnlasa/v112y2017i517p446-456.html

## individual test items
get_model_item <- function(outcome_item){
    return(lmer(as.formula(gsub("gcp", outcome_item, m0_form)), data = dt))
}

m20 <- get_model_item("memory"); m21 <- get_model_item("executive"); m22 <- get_model_item("language"); m23 <- get_model_item("visuospatial")
m24 <- get_model_item("wr_delayed"); m25 <- get_model_item("lm_delayed"); m26 <- get_model_item("animals"); m27 <- get_model_item("con_praxis"); m28 <- get_model_item("ravens")

## multiple imputation sensitivity model
mi_dt <- read_rds(required_input_file(file.path(derived_dir, "mi_data.rds")))
format_mi <- function(imp_num){
    data <- as.data.table(complete(mi_dt, imp_num))
    data <- melt.data.table(data, id.vars = c("prim_key", "educ", model_covs[!model_covs == "age"]), variable.name = "wave", measure.vars = patterns("gcp", "age", "time"), value.name = c("gcp", "age", "time"))
    return(as.data.frame(data))
}
mi_dt <- miceadds::datalist2mids(lapply(1:mi_dt$m, format_mi))
m29 <- pool(with(mi_dt, lmer(as.formula(m0_form))))

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

model_results <- rbindlist(lapply(0:29, get_coefficients))

# FORMAT GRAPH DATA -----------------------------------------------------------

## make parameter names consistent
model_results[, parameter := gsub("_L1$", "", parameter)][, parameter := gsub("\\.", " ", parameter)][parameter == "educMiddlesecondary school:time", parameter := "educMiddle-secondary school:time"] 
model_results[, parameter := fcase(parameter == "educLessthanprimary:time", "educLess than primary:time", 
                                   parameter == "educPrimaryschool:time", "educPrimary school:time", 
                                   parameter == "educMiddlesecondaryschool:time", "educMiddle-secondary school:time", 
                                   parameter == "educHighersecondaryschoolandup:time", "educHigher secondary school and up:time", 
                                   default = as.character(parameter))]
## get education categories
model_results[, educ_cat := gsub("educ", "", parameter)][, educ_cat := gsub(":time", "", educ_cat)][, educ_cat := factor(educ_cat, levels = c("Less than primary", "Primary school", "Middle-secondary school", "Higher secondary school and up"))]

## merge on labels and categories
label_dt <- data.table(model_num = model_results[, unique(model_num)], 
                       labels = c("Mixed effects model (base)", "Weighted GEE model (death)", "Weighted GEE model (survey + attrition)", 
                                  "Weighted GEE model (survey + death + attrition)", "Joint model", "Age (60-69)", "Age (70-79)", 
                                  "Age (80+)", "Urbanicity (rural)", "Urbanicity (urban)", "Gender (men)", "Gender (women)",
                                  "W1 cognition (below median W1 cog)", "W1 cognition (above median W1 cog)", 
                                  "W2 follow-up (below median follow-up)", "W2 follow-up (above median follow-up)", 
                                  "W1 CDR (normal)", "W1 CDR (MCI/Questionable)", "W1 CDR (dementia)",
                                  "10th percentile", "25th percentile", "50th percentile", "75th percentile", "90th percentile",
                                  "Memory", "Executive functioning", "Language", "Visuospatial functioning", 
                                  "Delayed word recall", "Delayed story recall (logical memory)", "Animal naming", "Constructional praxis", "Raven's progressive matrices", 
                                  "Mixed effects model \n(multiple imputation)"), 
                       category = c(rep("Base",5), rep("Stratified", 14), rep("Quantiles", 5), rep("Other cognitive outcomes", 9), "Base"))
label_dt[, labels := factor(labels, levels = rev(label_dt[, labels]))]                       
model_results <- merge(model_results, label_dt, by = "model_num", sort = FALSE)

## results section 
model_results[labels == "Mixed effects model (base)"]

# MAKE GRAPHS -----------------------------------------------------------

## initial sensitivities
graph1_models <- c(0,5:11,paste0("19 - ", c(0.1,0.25,0.5,0.75,0.9)), 20:28) ## models to include in the first graph 

get_graphs <- function(c){
    y_limits <- c(-0.12, 0.01)
    dt_plot <- copy(model_results[category == c & model_num %in% graph1_models])
    dt_plot[, lower_cap := pmax(lower, y_limits[1])]
    dt_plot[, upper_cap := pmin(upper, y_limits[2])]
    dt_plot[, label_text := ifelse(lower_cap == y_limits[1] | upper_cap == y_limits[2], sprintf("%.2f (%.2f, %.2f)", estimate, lower, upper), "")]
    dt_plot[, label_text := gsub("\\-0\\.00", "0.00", label_text)]
    dt_plot[, yval := ifelse(lower_cap == y_limits[1], lower_cap + 0.03, upper_cap - 0.03)]
    

    plot <- ggplot(model_results[category == c & model_num %in% graph1_models], aes(x = labels, y = estimate, ymin = lower, ymax = upper, color = educ_cat)) +
        geom_point() +
        geom_segment(data = dt_plot[lower < y_limits[1]], aes(x = labels, xend = labels, y = upper_cap, yend = lower_cap, color = educ_cat),
                     arrow = arrow(length = unit(0.2, "cm"), type = "closed"), show.legend = FALSE) +
        geom_segment(data = dt_plot[upper > y_limits[2]], aes(x = labels, xend = labels, y = lower_cap, yend = upper_cap, color = educ_cat),
                     arrow = arrow(length = unit(0.2, "cm"), type = "closed"), show.legend = FALSE) +      
        #geom_text(data = dt_plot, aes(x = labels, y = yval, label = label_text), size = 3.5, position = position_nudge(x = 0.3), color = "black") +                                  
        geom_errorbar(width = 0.1) + 
        geom_hline(yintercept = 0, linetype = "dashed") +
        facet_wrap(~educ_cat, nrow = 1) + 
        coord_flip() + 
        labs(x = "", y = "") +
        scale_y_continuous(limits = c(-0.12,0.01), breaks = seq(-.1,0,0.025), oob = oob_keep, expand = expansion(0)) +
        scale_color_manual(name = "", values = c("#5CB85C", "#357EBD", "#D43F3A", "#9632B8")) + 
        theme_bw() +
        theme(legend.position = "none", plot.margin = unit(c(0,5.5,0,5.5), "pt"))
    if (c == "Other cognitive outcomes"){
        plot <- plot + labs(y = "Difference in rate of annual cognitive decline (Ref: No schooling)") + theme(legend.position = "bottom")
    }
    if (!c == "Base"){
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
    plot_layout(nrow = 4, widths = c(1, 0.1), heights = c(0.2,1.15,0.9,1.5))

ggsave(file.path(plot_dir, paste0("regression_comparison_", date, ".pdf")), plot= full_plot, height = 5, width = 12)

## show all estimates in a table 
model_table <- copy(model_results[model_num %in% graph1_models])
model_table[, cell := sprintf("%.3f (%.3f, %.3f)", estimate, lower, upper)]
model_table <- dcast(model_table, category + labels ~ educ_cat, value.var = "cell")
model_table <- model_table[order(-labels)]
write.xlsx(model_table, file.path(appendix_dir, paste0("model_estimates_table_", date, ".xlsx")))

## mortality models for appendix
mort_results <- copy(model_results[category == "Base" & model_num %in% c(0,1,4)])
mort_results[, labels := forcats::fct_recode(labels, 
                                              "Mixed effects\nmodel (base)" = "Mixed effects model (base)", 
                                              "Weighted GEE\nmodel (death)" = "Weighted GEE model (death)")]

mortplot <- ggplot(mort_results, aes(x = labels, y = estimate, ymin = lower, ymax = upper, color = educ_cat)) +
        geom_point() +
        geom_errorbar(width = 0.1) + 
        geom_hline(yintercept = 0, linetype = "dashed") +
        facet_wrap(~educ_cat, nrow = 1) + 
        coord_flip() + 
        labs(x = "", y = "Difference in rate of annual cognitive decline (Ref: No schooling)") +
        scale_y_continuous(limits = c(-0.102,0), breaks = seq(-.1,0,0.025), oob = oob_keep) +
        scale_color_manual(name = "", values = c("#5CB85C", "#357EBD", "#D43F3A", "#9632B8")) + 
        theme_bw() +
        theme(legend.position = "bottom", plot.margin = unit(c(0,5.5,0,5.5), "pt"), axis.text.x = element_text(size = 8))

ggsave(file.path(appendix_dir, paste0("mortality_regression_comparison_", date, ".pdf")), plot= mortplot, height = 3, width = 10)

## regression to the mean models for appendix 
regmean <- ggplot(model_results[model_num %in% as.character(c(12:13,16:18))], aes(x = labels, y = estimate, ymin = lower, ymax = upper, color = educ_cat)) +
        geom_point() +
        geom_errorbar(width = 0.1) + 
        geom_hline(yintercept = 0, linetype = "dashed") +
        facet_wrap(~educ_cat, nrow = 1) + 
        coord_flip() + 
        labs(x = "", y = "Difference in rate of annual cognitive decline (Ref: No schooling)") +
        scale_y_continuous(limits = c(-0.102,0), breaks = seq(-.1,0,0.025), oob = oob_keep) +
        scale_color_manual(name = "", values = c("#5CB85C", "#357EBD", "#D43F3A", "#9632B8")) + 
        theme_bw() +
        theme(legend.position = "bottom", plot.margin = unit(c(0,5.5,0,5.5), "pt"))

ggsave(file.path(appendix_dir, paste0("regmean_regression_comparison_", date, ".pdf")), plot= regmean, height = 3, width = 14)

## multiple imputation for appendix
mi_plot <- ggplot(model_results[model_num %in% as.character(c(0,29))], aes(x = labels, y = estimate, ymin = lower, ymax = upper, color = educ_cat)) +
        geom_point() +
        geom_errorbar(width = 0.1) + 
        geom_hline(yintercept = 0, linetype = "dashed") +
        facet_wrap(~educ_cat, nrow = 1) + 
        coord_flip() + 
        labs(x = "", y = "Difference in rate of annual cognitive decline (Ref: No schooling)") +
        scale_y_continuous(limits = c(-0.102,0), breaks = seq(-.1,0,0.025), oob = oob_keep) +
        scale_color_manual(name = "", values = c("#5CB85C", "#357EBD", "#D43F3A", "#9632B8")) + 
        theme_bw() +
        theme(legend.position = "bottom")

ggsave(file.path(appendix_dir, paste0("mi_regression_comparison_", date, ".pdf")), plot= mi_plot, height = 3, width = 14)


# EXPLORE CHANGE MODELS -----------------------------------------------------------

change_dt <- copy(dt) 

## data manipulations
change_dt[, N := .N, by = "prim_key"]
change_dt[wave == 1, w1_age := age]
change_dt[, w1_age := mean(w1_age, na.rm = TRUE), by = "prim_key"]

## change covs to be w1 age 
change_covs <- model_covs
change_covs <- gsub("age", "w1_age", change_covs)

## create wide data 
change_dt <- dcast(as.formula(paste0("prim_key + educ + ", paste0(change_covs, collapse = " + "), " ~ wave")), value.var = "gcp",  data = change_dt[N==2])

## create change variable 
change_dt[, change := `2` - `1`]

## create baseline minus baseline mean 
change_dt[, baseline_demeaned := `1` - mean(`1`, na.rm = TRUE)]

## run models 
change_model <- lm(paste0("change ~  educ + ", paste0(change_covs, collapse = " + ")), data = change_dt)
change_model_adj <- lm(paste0("change ~ `1` + educ + ", paste0(change_covs, collapse = " + ")), data = change_dt)
change_model_demeaned <- lm(paste0("change ~ baseline_demeaned + educ + ", paste0(change_covs, collapse = " + ")), data = change_dt) ## this is just the same as adjusting for baseline???

secondtime_model_adj <- lm(paste0("`2` ~ `1` + educ + ", paste0(change_covs, collapse = " + ")), data = change_dt)
