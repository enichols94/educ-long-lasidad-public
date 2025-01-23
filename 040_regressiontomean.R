##########################################################################
### Author: Emma Nichols
### Date: 01/17/2025
### Project: LASIDAD Educationa and longitudinal change
### Purpose: Descriptive analyses - regression to the mean
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

## calculate change scores
plot_dt <- dcast(data = plot_dt, formula = prim_key + educ ~ wave, value.var = "gcp")   
plot_dt[, change := `2` - `1`]

## create means data
mean_dt <- plot_dt[, .(mean_baseline = mean(`1`, na.rm = TRUE), mean_change = mean(change, na.rm = TRUE)), by = "educ"]

## regression to the mean plot
rm_plot <- ggplot(plot_dt, aes(x = `1`, y = change, color = educ, fill = educ)) + 
    geom_point(alpha = 0.4) + 
    geom_abline(intercept = c(-1,0,1), slope = 1, linetype = "dashed") +
    geom_point(data = mean_dt, aes(x = mean_baseline, y = mean_change), color = "black", shape = 23, size = 3) +
    scale_color_manual(name = "", values = c("#EEA236", "#5CB85C", "#357EBD", "#D43F3A", "#9632B8")) +
    scale_fill_manual(name = "", values = c("#EEA236", "#5CB85C", "#357EBD", "#D43F3A", "#9632B8")) +
    labs(x = "Baseline cognitive functioning", y = "Change in cognitive functioning") +
    theme_bw() + 
    theme(legend.position = "bottom")

ggsave(paste0(plot_dir, "regression_tomean_", date, ".pdf"), plot = rm_plot, width = 8, height = 6)

