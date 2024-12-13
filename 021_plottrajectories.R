##########################################################################
### Author: Emma Nichols
### Date: 12/12/2024
### Project: LASIDAD Educationa and longitudinal change
### Purpose: Descriptive analyses - plot trajectories
##########################################################################

rm(list = ls())

pacman::p_load(data.table, openxlsx, haven, readr, dplyr, magrittr, stringr,
               labelled, scales, gridExtra, grid, ggplot2, lubridate, JM, 
               survival, patchwork, survey, srvyr)
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

data <- read_rds(paste0(derived_dir, "processed_data_weights.rds"))
survival_dt <- data$survival
longitudinal_dt <- data$longitudinal

# PLOTTING --------------------------------------------------------------------

## keep people with two visits
longitudinal_dt[, num_visits := .N, by = "prim_key"]
plot_dt <- copy(longitudinal_dt[num_visits == 2])

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

ggsave(paste0(plot_dir, "descriptive_trajectories_", date, ".pdf"), plot = traj_plot_facet, width = 12, height = 6)

# PROPORTION IMPROVING ---------------------------------------------------------

improve_dt <- copy(plot_dt)
improve_dt <- dcast(data = improve_dt, formula = prim_key + educ ~ wave, value.var = "gcp")
improve_dt[, difference := `2` - `1`]

improve_dt[, improved := as.numeric(difference > 0)]
improve_dt[, mean(improved, na.rm = TRUE), by = "educ"][order(educ)]
