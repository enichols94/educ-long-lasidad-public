# Plot observed cognitive trajectories by education category.
#
# Input: data/derived/processed_data_weights.rds from 011_weights.R.
# Output: a date-stamped descriptive trajectory PDF.

source(file.path("R", "project_paths.R"))
load_required_packages(c(
  "data.table", "openxlsx", "haven", "readr", "dplyr", "magrittr",
  "stringr", "labelled", "scales", "gridExtra", "grid", "ggplot2",
  "lubridate", "JM", "survival", "patchwork", "survey", "srvyr",
  "forcats"
))
date <- gsub("-", "_", Sys.Date())
set.seed(6541)

# PATHS -------------------------------------------------------------------

paths <- project_paths()
derived_dir <- paths$derived
plot_dir <- paths$trajectory_figures
ensure_output_directory(plot_dir)

# READ DATA ------------------------------------------------------------------

data <- read_rds(required_input_file(file.path(derived_dir, "processed_data_weights.rds")))
survival_dt <- data$survival
longitudinal_dt <- data$longitudinal

# RESTRICTIONS ---------------------------------------------------------------

## number lost to mortality and loss to follow-up  
survival_dt[, sum(died)]; survival_dt[, sum(as.numeric(refused == 1 | attrited == 1))]

## remove missing data
model_covs <- c("age", "gender", "caste", "rural", "educ_dad", "childhood_finance", "childhood_health")
dt <- copy(longitudinal_dt); surv_dt <- copy(survival_dt)
for (var in c("gcp", "educ", model_covs)){
    message(paste0("Missing individuals in ", var, ": ", dt[is.na(get(var)), length(unique(prim_key))], " (", 
                   dt[is.na(get(var)) & wave == 2, length(unique(prim_key))], " longitudinal)"))
    dt <- dt[!is.na(get(var))]
    message(paste0("Now excluding N=", nrow(surv_dt[!prim_key %in% dt[, unique(prim_key)]]), " from data completely"))
}
surv_dt <- surv_dt[prim_key %in% dt[, unique(prim_key)]] ## only keep those in longitudinal data

# PLOTTING --------------------------------------------------------------------

## keep people with two visits with cognition
longitudinal_dt[!is.na(gcp), num_visits := .N, by = "prim_key"]
plot_dt <- copy(longitudinal_dt[num_visits == 2])
plot_dt[, educ := forcats::fct_recode(educ,
                                     "Higher secondary\nschool and up" = "Higher secondary school and up")]

traj_plot <- ggplot(plot_dt, aes(x = time, y = gcp, color = educ)) +
    geom_point(alpha = 0.1) +
    geom_line(alpha = 0.1, aes(group = prim_key)) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 3) +
    scale_color_manual(name = "", values = c("#EEA236", "#5CB85C", "#357EBD", "#D43F3A", "#9632B8")) +
    labs(x = "Time since baseline (yrs)", y = "General cognitive functioning") +
    theme_bw()

traj_plot_facet <- ggplot(plot_dt, aes(x = time, y = gcp, color = educ)) +
    geom_point(alpha = 0.1) +
    geom_line(alpha = 0.1, aes(group = prim_key)) +
    facet_wrap(~educ, nrow = 1) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 3) +
    scale_color_manual(name = "", values = c("#EEA236", "#5CB85C", "#357EBD", "#D43F3A", "#9632B8")) +
    labs(x = "Time since baseline (yrs)", y = "General cognitive functioning") +
    theme_bw() +
    theme(legend.position = "bottom")

ggsave(file.path(plot_dir, paste0("descriptive_trajectories_", date, ".pdf")), plot = traj_plot_facet, width = 9, height = 4)

# PROPORTION IMPROVING ---------------------------------------------------------

improve_dt <- copy(plot_dt)
improve_dt <- dcast(data = improve_dt, formula = prim_key + educ ~ wave, value.var = "gcp")
improve_dt[, difference := `2` - `1`]

improve_dt[, improved := as.numeric(difference > 0)]
improve_dt[, mean(improved, na.rm = TRUE), by = "educ"][order(educ)]

# COGNITIVE LEVEL AT WAVE 2 ----------------------------------------------------

plot_dt[wave == 2, mean(gcp), by = "educ"][order(educ)]

