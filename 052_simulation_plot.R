# Assemble the final four-panel simulation figure.
#
# Inputs: simulation1.rds from script 050, simulation2.rds from script 051,
# and the two user-provided Educ-cogdecline simulation JPEG panels.
# Output: a date-stamped simulation PDF.

source(file.path("R", "project_paths.R"))
load_required_packages(c(
  "data.table", "openxlsx", "haven", "readr", "dplyr", "magrittr",
  "stringr", "labelled", "scales", "gridExtra", "grid", "ggplot2",
  "lubridate", "JM", "figpatch", "survival", "patchwork", "survey",
  "srvyr", "doBy", "lme4", "pbapply", "future.apply"
))
date <- gsub("-", "_", Sys.Date())
set.seed(6541)

# PATHS -------------------------------------------------------------------

paths <- project_paths()
rawdata_dir <- paths$source
derived_dir <- paths$derived
plot_dir <- paths$simulation_figures
ensure_output_directory(plot_dir)

# READ FILES --------------------------------------------------------------

sim1_fig <- read_rds(required_input_file(file.path(derived_dir, "simulation1.rds")))
sim2_fig <- read_rds(required_input_file(file.path(derived_dir, "simulation2.rds")))

sim1_structure <- figpatch::fig(required_input_file(file.path(rawdata_dir, "Educ-cogdecline-sim1.jpg")))
sim2_structure <- figpatch::fig(required_input_file(file.path(rawdata_dir, "Educ-cogdecline-sim2.jpg")))

# MAKE PLOT ----------------------------------------------------------------

simplot <- sim1_structure + sim1_fig + sim2_structure + sim2_fig +
    plot_layout(ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1)) + 
    plot_annotation(tag_levels = 'A', tag_suffix = ")") 

ggsave(file.path(plot_dir, paste0("simulation_plot_", date, ".pdf")), plot = simplot, width = 7, height = 5)

