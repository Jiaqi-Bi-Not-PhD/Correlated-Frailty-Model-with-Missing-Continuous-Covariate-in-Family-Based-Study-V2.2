## Author: Jiaqi Bi

## Run the paper-facing continuous-covariate simulation grid directly from R.
## Each scenario is launched as a separate R process, and the corresponding
## summary script runs after all scenarios finish.

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_file <- if (length(file_arg)) {
  sub("^--file=", "", file_arg[1])
} else {
  "scripts/run_simulation_grid.R"
}
script_file <- normalizePath(script_file, mustWork = TRUE)
repository_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)

parse_options <- function(args) {
  out <- list()
  for (arg in args) {
    if (!grepl("^--[^=]+=", arg)) {
      stop("Arguments must use --name=value syntax; got: ", arg, call. = FALSE)
    }
    key <- sub("^--([^=]+)=.*$", "\\1", arg)
    key <- gsub("-", "_", key, fixed = TRUE)
    out[[key]] <- sub("^--[^=]+=", "", arg)
  }
  out
}

options <- parse_options(commandArgs(trailingOnly = TRUE))
`%||%` <- function(x, y) if (is.null(x)) y else x
value_or <- function(name, default) options[[name]] %||% default

component <- tolower(value_or("component", "all"))
valid_components <- c("all", "no-missing", "comparators", "pdmi")
if (!component %in% valid_components) {
  stop(
    "--component must be one of: ",
    paste(valid_components, collapse = ", "),
    call. = FALSE
  )
}

as_positive_integer <- function(value, name) {
  out <- suppressWarnings(as.integer(value))
  if (is.na(out) || out < 1L) {
    stop(name, " must be a positive integer; got: ", value, call. = FALSE)
  }
  out
}

n_replicates <- as_positive_integer(value_or("n", "1000"), "--n")
n_families <- as_positive_integer(value_or("families", "498"), "--families")
n_cores <- as_positive_integer(value_or("cores", "1"), "--cores")
m_pdmi <- as_positive_integer(value_or("m_pdmi", "20"), "--m-pdmi")
pdmi_numit <- as_positive_integer(value_or("pdmi_numit", "10"), "--pdmi-numit")
m_smcfcs <- as_positive_integer(value_or("m_smcfcs", "20"), "--m-smcfcs")
run_prefix <- value_or("run_prefix", "reviewer")
force <- tolower(value_or("force", "false")) %in% c("1", "true", "yes", "y")
dry_run <- tolower(value_or("dry_run", "false")) %in% c("1", "true", "yes", "y")
output_root <- normalizePath(
  value_or("output", file.path(repository_root, "reproduced-results")),
  mustWork = FALSE
)
if (!dry_run) dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

maxit <- options$maxit
if (!is.null(maxit)) maxit <- as_positive_integer(maxit, "--maxit")

restore_environment <- function(old_values) {
  missing_before <- is.na(old_values)
  if (any(missing_before)) Sys.unsetenv(names(old_values)[missing_before])
  present_before <- !missing_before
  if (any(present_before)) {
    do.call(Sys.setenv, as.list(old_values[present_before]))
  }
}

run_r_script <- function(relative_script, env) {
  script <- file.path(repository_root, relative_script)
  if (!file.exists(script)) stop("Missing R script: ", script, call. = FALSE)
  if (dry_run) {
    message(
      "[dry run] Rscript ", relative_script, " | ",
      paste(paste0(names(env), "=", env), collapse = " ")
    )
    return(invisible(0L))
  }
  keys <- names(env)
  old_values <- Sys.getenv(keys, unset = NA_character_)
  names(old_values) <- keys
  on.exit(restore_environment(old_values), add = TRUE)
  env_values <- as.list(as.character(env))
  names(env_values) <- keys
  do.call(Sys.setenv, env_values)
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    args = shQuote(script)
  )
  if (!identical(status, 0L)) {
    stop("Simulation script failed: ", relative_script, call. = FALSE)
  }
  invisible(status)
}

base_env <- c(
  SIM_N = as.character(n_replicates),
  SIM_NFAM = as.character(n_families),
  SIM_REP_START = "1",
  SIM_REP_END = as.character(n_replicates),
  SIM_CORES = as.character(n_cores)
)
if (!is.null(maxit)) base_env <- c(base_env, SIM_MAXIT = as.character(maxit))

run_no_missing <- function() {
  results_root <- file.path(output_root, "no-missing")
  run_label <- paste0(run_prefix, "_no_missing_n", n_replicates)
  for (sigma_u2 in c(0.2, 0.5)) {
    message("No-missing: sigma_u2=", sigma_u2)
    run_r_script(
      "code/comparators-may23/No Missing/run_no_missing_benchmark.R",
      c(
        base_env,
        SIM_SIGMA_U2 = as.character(sigma_u2),
        SIM_RESULTS_ROOT = results_root,
        SIM_RUN_LABEL = run_label,
        SIM_FORCE_BENCHMARKS = if (force) "1" else "0"
      )
    )
  }
  run_r_script(
    "code/comparators-may23/Misc/summarize_results.R",
    c(SIM_RESULTS_ROOT = results_root, SIM_RUN_LABEL = run_label)
  )
}

run_comparators <- function() {
  results_root <- file.path(output_root, "comparators")
  run_label <- paste0(run_prefix, "_comparators_n", n_replicates, "_m", m_smcfcs)
  for (sigma_u2 in c(0.2, 0.5)) {
    for (missing_rate in c(0.20, 0.50, 0.80)) {
      for (method in c("cca", "smcfcs")) {
        message(
          "Comparator: method=", method,
          " sigma_u2=", sigma_u2,
          " missing=", missing_rate
        )
        run_r_script(
          "code/comparators-may23/Continuous Only/run_continuous_simulation.R",
          c(
            base_env,
            SIM_SIGMA_U2 = as.character(sigma_u2),
            SIM_MISS = as.character(missing_rate),
            SIM_METHOD = method,
            SIM_M_SMCFCS = as.character(m_smcfcs),
            SIM_RESULTS_ROOT = results_root,
            SIM_RUN_LABEL = run_label,
            SIM_FORCE_RESULTS = if (force) "1" else "0"
          )
        )
      }
    }
  }
  run_r_script(
    "code/comparators-may23/Misc/summarize_results.R",
    c(SIM_RESULTS_ROOT = results_root, SIM_RUN_LABEL = run_label)
  )
}

run_pdmi <- function() {
  results_root <- file.path(output_root, "pdmi")
  run_label <- paste0(run_prefix, "_pdmi_n", n_replicates, "_m", m_pdmi)
  for (sigma_u2 in c(0.2, 0.5)) {
    for (missing_rate in c(0.20, 0.50, 0.80)) {
      for (method in c("c-o-pdmi", "c-r-pdmi")) {
        message(
          "PDMI: method=", method,
          " sigma_u2=", sigma_u2,
          " missing=", missing_rate
        )
        run_r_script(
          "code/pdmi-v2.2/Continuous Only/run_continuous_pdmi_simulation.R",
          c(
            base_env,
            SIM_SIGMA_U2 = as.character(sigma_u2),
            SIM_MISS = as.character(missing_rate),
            SIM_METHOD = method,
            SIM_M = as.character(m_pdmi),
            SIM_PDMI_NUMIT = as.character(pdmi_numit),
            SIM_RESULTS_ROOT = results_root,
            SIM_RUN_LABEL = run_label,
            SIM_FORCE_RESULTS = if (force) "1" else "0"
          )
        )
      }
    }
  }
  run_r_script(
    "code/pdmi-v2.2/Misc/summarize_results.R",
    c(SIM_RESULTS_ROOT = results_root, SIM_RUN_LABEL = run_label)
  )
}

message(
  "Starting component=", component,
  " n=", n_replicates,
  " families=", n_families,
  " cores=", n_cores,
  " output=", output_root
)
if (component %in% c("all", "no-missing")) run_no_missing()
if (component %in% c("all", "comparators")) run_comparators()
if (component %in% c("all", "pdmi")) run_pdmi()
message("Requested simulation component completed.")
