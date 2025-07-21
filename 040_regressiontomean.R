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
plot_dir <- paste0(dir, "paper/regmean_fig/")

# READ DATA ------------------------------------------------------------------

data <- read_rds(paste0(derived_dir, "processed_data_weights.rds"))
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

## keep people with two visits
longitudinal_dt[, num_visits := .N, by = "prim_key"]
plot_dt <- copy(longitudinal_dt[num_visits == 2])

## calculate change scores
plot_dt <- dcast(data = plot_dt, formula = prim_key + educ ~ wave, value.var = "gcp")   
plot_dt[, change := `2` - `1`]
plot_dt[, cor.test(`1`, change, use = "pairwise.complete.obs")] # correlation between baseline and change scores

## calculate the mean change and baseline value for those declining or improving more than 0.5 
plot_dt[, `:=` (improve = as.numeric(change > 0.5), decline = as.numeric(change < -0.5))]

## calculate the proportion in each education group for those who improve vs. decline 
improve_dt <- plot_dt[improve == 1, .N, by = educ][, prop_improve := N/sum(N)]
decline_dt <- plot_dt[improve == 0, .N, by = educ][, prop_decline := N/sum(N)]
change_dt <- merge(improve_dt[, .(educ, prop_improve)], decline_dt[, .(educ, prop_decline)], by = "educ", all = TRUE) # merge to get proportions of improve and decline
change_dt[, ratio := prop_decline/prop_improve]; print(change_dt)

## create decline data
decline_dt <- data.table(type = c("improve", "decline"), 
                         `1` = c(plot_dt[improve == 1, mean(`1`)], plot_dt[decline == 1, mean(`1`)]), # baseline cognitive functioning for those improving or declining
                         `2` = c(plot_dt[improve == 1, mean(`2`)], plot_dt[decline == 1, mean(`2`)])) # wave 2 cognitive functioning for those improving or declining
decline_dt

## create means data
mean_dt <- plot_dt[, .(mean_baseline = mean(`1`, na.rm = TRUE), mean_time2 = mean(`2`, na.rm = TRUE), mean_change = mean(change, na.rm = TRUE)), by = "educ"]
mean_dt

## regression to the mean plot 1
rm_plot1 <- ggplot(plot_dt, aes(x = `1`, y = `2`)) + 
    geom_point(alpha = 0.4, aes(color = educ, fill = educ)) + 
    geom_point(data = decline_dt, aes(x = `1`, y = `2`), shape = 13, size = 5) + ## highlight the decline/improve points
    geom_vline(xintercept = decline_dt[, `1`], linetype = "dashed", color = "black", alpha = 0.5) + ## vertical line at the baseline cognitive functioning of those declining or improving
    geom_abline(intercept = c(-0.5,0,0.5), slope = 1, linetype = "dotdash") + ## 0.9 is the standard deviation of gcp at time 1
    scale_color_manual(name = "", values = c("#EEA236", "#5CB85C", "#357EBD", "#D43F3A", "#9632B8")) +
    scale_fill_manual(name = "", values = c("#EEA236", "#5CB85C", "#357EBD", "#D43F3A", "#9632B8")) +
    labs(x = "Baseline cognitive functioning", y = "Wave 2 cognitive functioning") +
    guides(alpha = "none", fill = "none") +
    theme_bw() + 
    theme(legend.position = "bottom")

rm_plot2 <- ggplot(plot_dt, aes(x = `1`, y = change)) + 
    geom_point(alpha = 0.4, aes(color = educ, fill = educ)) + 
    geom_hline(yintercept = 0, linetype = "dotdash", color = "black") + 
    geom_smooth(method = "lm", se = FALSE, color = "black", aes(alpha = 0.2), linewidth = 0.5) + 
    geom_point(data = mean_dt, aes(x = mean_baseline, y = mean_change, color = educ, fill = educ), color = "black", shape = 23, size = 3) +
    scale_color_manual(name = "", values = c("#EEA236", "#5CB85C", "#357EBD", "#D43F3A", "#9632B8")) +
    scale_fill_manual(name = "", values = c("#EEA236", "#5CB85C", "#357EBD", "#D43F3A", "#9632B8")) +
    labs(x = "Wave 1 cognitive functioning", y = "Change in cognitive functioning") +
    guides(alpha = "none", fill = "none") +
    theme_bw() + 
    theme(legend.position = "bottom")

## black and white version 
rm_plot2_bw <- ggplot(plot_dt, aes(x = `1`, y = change)) + 
    geom_point(alpha = 0.4) + 
    geom_hline(yintercept = 0, linetype = "dotdash", color = "black") + 
    geom_smooth(method = "lm", se = FALSE, color = "black", aes(alpha = 0.2), linewidth = 0.5) + 
    labs(x = "Baseline cognitive functioning", y = "Change in cognitive functioning") +
    guides(alpha = "none", fill = "none") +
    theme_bw() + 
    theme(legend.position = "bottom")

## pull in simulation plot 
sim_plot <- read_rds(paste0(derived_dir, "simulated_rtm_plot.rds"))

rm_fullplot <- rm_plot2 + sim_plot +
    plot_annotation(tag_levels = ("A"), tag_suffix = ".") +
    plot_layout(nrow = 1, guides = "collect") & theme(legend.position = "bottom")     

ggsave(paste0(plot_dir, "regression_tomean_", date, ".pdf"), plot = rm_fullplot, width = 14, height = 6)

#ggsave(paste0(plot_dir, "regression_tomean_bw_melodem_", date, ".pdf"), plot = rm_plot2_bw, width = 7, height = 6)
