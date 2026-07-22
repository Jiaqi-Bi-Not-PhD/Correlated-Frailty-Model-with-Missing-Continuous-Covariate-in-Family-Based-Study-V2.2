## Author: Jiaqi Bi

## ============================================================
## No-missing full-data frailtypack benchmark.
## Environment variables:
##   SIM_N, SIM_NFAM, SIM_SIGMA_U2, SIM_REP_START, SIM_REP_END,
##   SIM_RESULTS_ROOT, SIM_RUN_LABEL
## ============================================================

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "No Missing/run_no_missing_benchmark.R"
script_file <- gsub("~+~", " ", script_file, fixed = TRUE)
code_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = FALSE)
if (!dir.exists(file.path(code_root, "Shared"))) code_root <- normalizePath(getwd(), mustWork = FALSE)
source(file.path(code_root, "Shared", "source_all.R"))

v2_config_from_env <- function() {
  cfg <- v2_default_config()
  cfg$n_families <- as.integer(Sys.getenv("SIM_NFAM", cfg$n_families))
  cfg$B_sim <- as.integer(Sys.getenv("SIM_N", cfg$B_sim))
  cfg$results_root <- Sys.getenv("SIM_RESULTS_ROOT", cfg$results_root)
  cfg$run_label <- Sys.getenv("SIM_RUN_LABEL", cfg$run_label)
  cfg$frailtypack_maxit <- as.integer(Sys.getenv("SIM_MAXIT", cfg$frailtypack_maxit))
  cfg$skip_existing_benchmarks <- !identical(Sys.getenv("SIM_FORCE_BENCHMARKS", "0"), "1")
  cfg
}

v2_run_no_missing_one <- function(replicate_id, sigma_u2, config) {
  seeds <- v2_seed_streams(replicate_id, "none", NA_real_, "full-data")
  metadata <- v2_result_metadata(replicate_id, "none", sigma_u2, NA_real_,
                                 "full-data", "none", 0L, seeds, config)
  existing <- v2_existing_benchmark_path(metadata, config)
  if (!is.na(existing)) return(existing)
  out <- tryCatch({
    gen <- v2_generate_complete_data(replicate_id, sigma_u2, seeds$complete_seed, config)
    fit <- v2_fit_frailtypack(gen$dat, gen$K, config,
                              init_omega = v2_actual_omega(sigma_u2, config),
                              method = "full-data")
    pen <- if (isTRUE(fit$convergence)) {
      v2_penetrance_from_fit(fit$omega, fit$vcov_omega, config)
    } else data.frame()
    v2_pack_result(
      metadata,
      fit = fit,
      penetrance = pen,
      diagnostics = list(generation_attempts = gen$generation_attempts,
                         n_analysis = nrow(gen$dat),
                         n_families = length(unique(gen$dat$famID))),
      convergence = fit$convergence,
      failure_reason = fit$failure_reason
    )
  }, error = function(e) {
    v2_pack_result(metadata, convergence = FALSE, failure_reason = conditionMessage(e))
  })
  v2_write_result(out, config)
}

config <- v2_config_from_env()
v2_require_packages(include_frailtypack = TRUE)
v2_load_family_sources(code_root)

sigma_u2 <- v2_parse_single_numeric_env("SIM_SIGMA_U2", config$sigma_u2_grid[1],
                                        "SIM_SIGMA_U2")
rep_start <- as.integer(Sys.getenv("SIM_REP_START", "1"))
rep_end <- as.integer(Sys.getenv("SIM_REP_END", as.character(config$B_sim)))
reps <- rep_start:rep_end
cores <- v2_get_job_cores()

message("V2 no-missing benchmark: reps=", length(reps), " sigma=", sigma_u2,
        " cores=", cores)
config$sigma_u2_grid <- sigma_u2
paths <- parallel::mclapply(reps, v2_run_no_missing_one, sigma_u2 = sigma_u2, config = config,
                            mc.cores = cores, mc.preschedule = TRUE, mc.silent = TRUE)
dir.create(file.path(config$results_root, config$run_label), recursive = TRUE, showWarnings = FALSE)
saveRDS(paths, file.path(config$results_root, config$run_label,
                         paste0("no_missing_paths_", v2_number_tag("sigma_u2", sigma_u2), ".rds")))
message("Completed no-missing benchmark.")
