# Shared project configuration for the public analysis scripts.
#
# Scripts are expected to be launched from the repository root. Restricted
# datasets remain outside the repository and are located through environment
# variables; repository inputs and generated outputs use relative paths.

load_required_packages <- function(packages) {
  missing_packages <- packages[!vapply(
    packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )]

  if (length(missing_packages) > 0) {
    stop(
      "Install the required R package(s) before running this script: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(lapply(packages, library, character.only = TRUE))
}

project_paths <- function() {
  project_root <- normalizePath(
    getwd(),
    winslash = "/",
    mustWork = TRUE
  )
  helper_file <- file.path(project_root, "R", "project_paths.R")

  if (!file.exists(helper_file)) {
    stop(
      "Run this script from the repository root (<PROJECT_ROOT>).",
      call. = FALSE
    )
  }

  list(
    root = project_root,
    source = file.path(project_root, "data", "source"),
    derived = file.path(project_root, "data", "derived"),
    plots = file.path(project_root, "plots"),
    appendix = file.path(project_root, "paper", "appendix_figs"),
    model_figures = file.path(project_root, "paper", "model_sensitivities_fig"),
    practice_figures = file.path(project_root, "paper", "practiceeffects_fig"),
    simulation_figures = file.path(project_root, "paper", "simulation_fig"),
    table1 = file.path(project_root, "paper", "table1"),
    trajectory_figures = file.path(project_root, "paper", "trajectories_fig")
  )
}

required_environment_directory <- function(variable_name) {
  directory <- Sys.getenv(variable_name, unset = "")

  if (!nzchar(directory)) {
    stop(
      "Set ", variable_name, " to the required controlled-data directory.",
      call. = FALSE
    )
  }
  if (!dir.exists(directory)) {
    stop(variable_name, " does not identify an existing directory.", call. = FALSE)
  }

  normalizePath(directory, winslash = "/", mustWork = TRUE)
}

required_input_file <- function(path) {
  if (!file.exists(path)) {
    stop("Required input file is missing: ", basename(path), call. = FALSE)
  }

  path
}

ensure_output_directory <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(path)) {
    stop("Unable to create the required output directory.", call. = FALSE)
  }

  invisible(path)
}
