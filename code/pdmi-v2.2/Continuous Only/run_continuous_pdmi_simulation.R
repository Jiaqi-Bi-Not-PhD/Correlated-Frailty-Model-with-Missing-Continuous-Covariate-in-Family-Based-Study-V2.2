## Author: Jiaqi Bi

## ============================================================
## Continuous-PRS missingness simulations for V2.2 PDMI.
## Methods: all, c-o-pdmi, c-r-pdmi. MI-SMCFCS is not rerun in V2.2.
## Environment variables:
##   SIM_N, SIM_NFAM, SIM_SIGMA_U2, SIM_MISS, SIM_METHOD,
##   SIM_REP_START, SIM_REP_END, SIM_RESULTS_ROOT, SIM_RUN_LABEL,
##   SIM_M, SIM_PDMI_NUMIT
## ============================================================

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "Continuous Only/run_continuous_pdmi_simulation.R"
script_file <- gsub("~+~", " ", script_file, fixed = TRUE)
code_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = FALSE)
if (!dir.exists(file.path(code_root, "Shared"))) code_root <- normalizePath(getwd(), mustWork = FALSE)
source(file.path(code_root, "Shared", "source_all.R"))

v22_config_from_env <- function() {
  cfg <- v22_default_config()
  cfg$n_families <- as.integer(Sys.getenv("SIM_NFAM", cfg$n_families))
  cfg$B_sim <- as.integer(Sys.getenv("SIM_N", cfg$B_sim))
  cfg$M_imp_pdmi <- as.integer(Sys.getenv("SIM_M", cfg$M_imp_pdmi))
  cfg$pdmi_numit <- as.integer(Sys.getenv("SIM_PDMI_NUMIT", cfg$pdmi_numit))
  cfg$results_root <- Sys.getenv("SIM_RESULTS_ROOT", cfg$results_root)
  cfg$run_label <- Sys.getenv("SIM_RUN_LABEL", cfg$run_label)
  cfg$frailtypack_maxit <- as.integer(Sys.getenv("SIM_MAXIT", cfg$frailtypack_maxit))
  cfg$skip_existing_results <- !identical(Sys.getenv("SIM_FORCE_RESULTS", "0"), "1")
  cfg$skip_existing_benchmarks <- FALSE
  cfg
}

v22_continuous_methods <- function(method_env) {
  method_env <- tolower(method_env)
  if (method_env == "all") return(c("c-o-pdmi", "c-r-pdmi"))
  methods <- trimws(unlist(strsplit(method_env, ","), use.names = FALSE))
  if ("smcfcs" %in% methods) {
    stop(
      "V2.2 does not run MI-SMCFCS; use code/comparators-may23 ",
      "with SIM_M_SMCFCS=20."
    )
  }
  methods[nzchar(methods)]
}

v22_continuous_method_metadata <- function(method, replicate_id, sigma_u2, miss_rate, config) {
  seeds <- v22_seed_streams(replicate_id, "continuous", miss_rate, method)
  method_l <- tolower(method)
  if (method_l %in% c("c-o-pdmi", "c-o", "co", "oracle")) {
    return(v22_result_metadata(replicate_id, "continuous", sigma_u2, miss_rate,
                              "C-O-PDMI", "C-O", config$M_imp_pdmi, seeds, config))
  }
  if (method_l %in% c("c-r-pdmi", "c-r", "cr", "realistic")) {
    return(v22_result_metadata(replicate_id, "continuous", sigma_u2, miss_rate,
                              "C-R-PDMI", "C-R", config$M_imp_pdmi, seeds, config))
  }
  v22_result_metadata(replicate_id, "continuous", sigma_u2, miss_rate,
                     method, NA_character_, NA_integer_, seeds, config)
}

v22_pack_continuous_method <- function(method, replicate_id, sigma_u2, miss_rate,
                                      dat_miss, K, mask_diag, config, seeds,
                                      pedigree_dat = NULL) {
  method_l <- tolower(method)
  if (method_l %in% c("c-o-pdmi", "c-o", "co", "oracle")) {
    metadata <- v22_continuous_method_metadata(method, replicate_id, sigma_u2, miss_rate, config)
    existing <- v22_existing_result_path(metadata, config)
    if (!is.na(existing)) return(existing)
    out <- v22_run_continuous_pdmi(dat_miss, K, "C-O", config, seeds$method_seed)
    return(v22_pack_result(metadata, pooled = out$pooled, penetrance = out$penetrance,
                          diagnostics = c(list(mask = mask_diag), out$diagnostics %||% list()),
                          convergence = out$convergence,
                          failure_reason = out$failure_reason))
  }
  if (method_l %in% c("c-r-pdmi", "c-r", "cr", "realistic")) {
    metadata <- v22_continuous_method_metadata(method, replicate_id, sigma_u2, miss_rate, config)
    existing <- v22_existing_result_path(metadata, config)
    if (!is.na(existing)) return(existing)
    out <- v22_run_continuous_pdmi(dat_miss, K, "C-R", config, seeds$method_seed)
    return(v22_pack_result(metadata, pooled = out$pooled, penetrance = out$penetrance,
                          diagnostics = c(list(mask = mask_diag), out$diagnostics %||% list()),
                          convergence = out$convergence,
                          failure_reason = out$failure_reason))
  }
  stop("Unknown continuous V2.2 PDMI method: ", method)
}

v22_run_continuous_one <- function(replicate_id, sigma_u2, miss_rate, methods, config) {
  cfg <- config
  cfg$sigma_u2_grid <- sigma_u2
  preflight <- lapply(methods, function(method) {
    md <- v22_continuous_method_metadata(method, replicate_id, sigma_u2, miss_rate, cfg)
    existing <- v22_existing_result_path(md, cfg)
    list(method = method, metadata = md, existing = existing)
  })
  if (all(vapply(preflight, function(x) !is.na(x$existing), logical(1)))) {
    return(lapply(preflight, `[[`, "existing"))
  }
  methods <- vapply(preflight[vapply(preflight, function(x) is.na(x$existing), logical(1))],
                    `[[`, character(1), "method")
  gen_seeds <- v22_seed_streams(replicate_id, "continuous", miss_rate, "common")
  gen <- v22_generate_complete_data(replicate_id, sigma_u2, gen_seeds$complete_seed, cfg)
  dat_miss <- v22_apply_continuous_missingness(gen$dat, miss_rate, gen_seeds$missing_mask_seed, cfg)
  mask_diag <- v22_missing_diagnostics(dat_miss, miss_rate)
  lapply(methods, function(method) {
    seeds <- v22_seed_streams(replicate_id, "continuous", miss_rate, method)
    res <- tryCatch(v22_pack_continuous_method(method, replicate_id, sigma_u2, miss_rate,
                                              dat_miss, gen$K, mask_diag, cfg, seeds,
                                              pedigree_dat = gen$pedigree_dat),
                    error = function(e) {
                      md <- v22_continuous_method_metadata(method, replicate_id, sigma_u2, miss_rate, cfg)
                      v22_pack_result(md, diagnostics = list(mask = mask_diag),
                                     convergence = FALSE,
                                     failure_reason = conditionMessage(e))
                    })
    if (is.character(res) && length(res) == 1L && file.exists(res)) res else v22_write_result(res, cfg)
  })
}

config <- v22_config_from_env()
methods <- v22_continuous_methods(Sys.getenv("SIM_METHOD", "all"))
v22_require_packages(include_frailtypack = TRUE, include_smcfcs = FALSE)
v22_load_family_sources(code_root)

sigma_u2 <- v22_parse_single_numeric_env("SIM_SIGMA_U2", config$sigma_u2_grid[1],
                                        "SIM_SIGMA_U2")
miss_rate <- v22_parse_single_numeric_env("SIM_MISS", config$missing_rates[1],
                                         "SIM_MISS")
rep_start <- as.integer(Sys.getenv("SIM_REP_START", "1"))
rep_end <- as.integer(Sys.getenv("SIM_REP_END", as.character(config$B_sim)))
reps <- rep_start:rep_end
cores <- v22_get_job_cores()

message("V2.2 continuous PDMI simulation: methods=", paste(methods, collapse = ","),
        " reps=", length(reps), " sigma=", sigma_u2, " miss=", miss_rate,
        " M_PDMI=", config$M_imp_pdmi,
        " numit=", config$pdmi_numit, " cores=", cores)
config$sigma_u2_grid <- sigma_u2
config$missing_rates <- miss_rate
paths <- parallel::mclapply(reps, function(replicate_id) {
  v22_run_continuous_one(replicate_id, sigma_u2, miss_rate, methods, config)
}, mc.cores = cores, mc.preschedule = FALSE, mc.silent = TRUE)
dir.create(file.path(config$results_root, config$run_label), recursive = TRUE, showWarnings = FALSE)
saveRDS(paths, file.path(config$results_root, config$run_label,
                         paste0("continuous_paths_",
                                v22_number_tag("sigma_u2", sigma_u2), "_",
                                v22_rate_tag(miss_rate), "_",
                                v22_clean_tag(paste(methods, collapse = "_")), "_",
                                rep_start, "_", rep_end, ".rds")))
message("Completed V2.2 continuous PDMI simulation.")
