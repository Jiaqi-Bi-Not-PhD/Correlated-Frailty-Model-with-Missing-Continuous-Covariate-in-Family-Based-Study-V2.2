## ============================================================
## Continuous-PRS missingness simulations.
## Author: Jiaqi Bi
## Methods: all, cca, smcfcs.
## Environment variables:
##   SIM_N, SIM_NFAM, SIM_SIGMA_U2, SIM_MISS, SIM_METHOD,
##   SIM_REP_START, SIM_REP_END, SIM_RESULTS_ROOT, SIM_RUN_LABEL
## ============================================================

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "Continuous Only/run_continuous_simulation.R"
script_file <- gsub("~+~", " ", script_file, fixed = TRUE)
code_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = FALSE)
if (!dir.exists(file.path(code_root, "Shared"))) code_root <- normalizePath(getwd(), mustWork = FALSE)
source(file.path(code_root, "Shared", "source_all.R"))

v2_config_from_env <- function() {
  cfg <- v2_default_config()
  cfg$n_families <- as.integer(Sys.getenv("SIM_NFAM", cfg$n_families))
  cfg$B_sim <- as.integer(Sys.getenv("SIM_N", cfg$B_sim))
  cfg$M_imp_smcfcs <- as.integer(Sys.getenv("SIM_M_SMCFCS", cfg$M_imp_smcfcs))
  cfg$results_root <- Sys.getenv("SIM_RESULTS_ROOT", cfg$results_root)
  cfg$run_label <- Sys.getenv("SIM_RUN_LABEL", cfg$run_label)
  cfg$frailtypack_maxit <- as.integer(Sys.getenv("SIM_MAXIT", cfg$frailtypack_maxit))
  cfg$skip_existing_results <- !identical(Sys.getenv("SIM_FORCE_RESULTS", "0"), "1")
  cfg$skip_existing_benchmarks <- !identical(Sys.getenv("SIM_FORCE_BENCHMARKS", "0"), "1")
  cfg
}

v2_continuous_methods <- function(method_env) {
  method_env <- tolower(method_env)
  if (method_env == "all") return(c("cca", "smcfcs"))
  methods <- trimws(strsplit(method_env, ",", fixed = TRUE)[[1]])
  invalid <- setdiff(methods, c("cca", "smcfcs"))
  if (length(invalid)) {
    stop(
      "This continuous-covariate comparator release supports only cca and smcfcs; got: ",
      paste(invalid, collapse = ", "),
      call. = FALSE
    )
  }
  methods
}

v2_continuous_method_metadata <- function(method, replicate_id, sigma_u2, miss_rate, config) {
  seeds <- v2_seed_streams(replicate_id, "continuous", miss_rate, method)
  method_l <- tolower(method)
  if (method_l == "cca") {
    return(v2_result_metadata(replicate_id, "continuous", sigma_u2, miss_rate,
                              "CCA", "none", 0L, seeds, config))
  }
  if (method_l == "smcfcs") {
    return(v2_result_metadata(replicate_id, "continuous", sigma_u2, miss_rate,
                              "MI-SMCFCS", "none", config$M_imp_smcfcs, seeds, config))
  }
  v2_result_metadata(replicate_id, "continuous", sigma_u2, miss_rate,
                     method, NA_character_, NA_integer_, seeds, config)
}

v2_pack_continuous_method <- function(method, replicate_id, sigma_u2, miss_rate,
                                      dat_miss, K, mask_diag, config, seeds, pedigree_dat = NULL) {
  method_l <- tolower(method)
  if (method_l == "cca") {
    metadata <- v2_result_metadata(replicate_id, "continuous", sigma_u2, miss_rate,
                                   "CCA", "none", 0L, seeds, config)
    existing <- v2_existing_result_path(metadata, config)
    if (!is.na(existing)) return(existing)
    cca <- v2_run_cca(dat_miss, K, config,
                      init_omega = v2_actual_omega(sigma_u2, config))
    pen <- if (isTRUE(cca$fit$convergence)) {
      v2_penetrance_from_fit(cca$fit$omega, cca$fit$vcov_omega, config)
    } else data.frame()
    return(v2_pack_result(metadata, fit = cca$fit, penetrance = pen,
                          diagnostics = list(mask = mask_diag, cca = cca$diagnostics),
                          convergence = cca$fit$convergence,
                          failure_reason = cca$fit$failure_reason))
  }
  if (method_l == "smcfcs") {
    metadata <- v2_result_metadata(replicate_id, "continuous", sigma_u2, miss_rate,
                                   "MI-SMCFCS", "none", config$M_imp_smcfcs, seeds, config)
    existing <- v2_existing_result_path(metadata, config)
    if (!is.na(existing)) return(existing)
    out <- v2_run_smcfcs_comparator(
      dat_miss,
      K,
      config = config,
      seed = seeds$method_seed
    )
    return(v2_pack_result(metadata, pooled = out$pooled, penetrance = out$penetrance,
                          diagnostics = c(list(mask = mask_diag), out$diagnostics %||% list()),
                          convergence = out$convergence,
                          failure_reason = out$failure_reason))
  }
  stop("Unknown continuous method: ", method)
}

v2_run_continuous_one <- function(replicate_id, sigma_u2, miss_rate, methods, config) {
  cfg <- config
  cfg$sigma_u2_grid <- sigma_u2
  preflight <- lapply(methods, function(method) {
    md <- v2_continuous_method_metadata(method, replicate_id, sigma_u2, miss_rate, cfg)
    existing <- v2_existing_result_path(md, cfg)
    list(method = method, metadata = md, existing = existing)
  })
  if (all(vapply(preflight, function(x) !is.na(x$existing), logical(1)))) {
    return(lapply(preflight, `[[`, "existing"))
  }
  methods <- vapply(preflight[vapply(preflight, function(x) is.na(x$existing), logical(1))],
                    `[[`, character(1), "method")
  gen_seeds <- v2_seed_streams(replicate_id, "continuous", miss_rate, "common")
  gen <- v2_generate_complete_data(replicate_id, sigma_u2, gen_seeds$complete_seed, cfg)
  dat_miss <- v2_apply_continuous_missingness(gen$dat, miss_rate, gen_seeds$missing_mask_seed, cfg)
  mask_diag <- v2_missing_diagnostics(dat_miss, miss_rate)
  lapply(methods, function(method) {
    seeds <- v2_seed_streams(replicate_id, "continuous", miss_rate, method)
    res <- tryCatch(v2_pack_continuous_method(method, replicate_id, sigma_u2, miss_rate,
                                              dat_miss, gen$K, mask_diag, cfg, seeds,
                                              pedigree_dat = gen$pedigree_dat),
                    error = function(e) {
                      md <- v2_result_metadata(replicate_id, "continuous", sigma_u2, miss_rate,
                                               method, NA_character_, NA_integer_, seeds, cfg)
                      v2_pack_result(md, diagnostics = list(mask = mask_diag),
                                     convergence = FALSE,
                                     failure_reason = conditionMessage(e))
                    })
    if (is.character(res) && length(res) == 1L && file.exists(res)) res else v2_write_result(res, cfg)
  })
}

config <- v2_config_from_env()
v2_require_packages(include_frailtypack = TRUE)
v2_load_family_sources(code_root)

methods <- v2_continuous_methods(Sys.getenv("SIM_METHOD", "all"))
sigma_u2 <- v2_parse_single_numeric_env("SIM_SIGMA_U2", config$sigma_u2_grid[1],
                                        "SIM_SIGMA_U2")
miss_rate <- v2_parse_single_numeric_env("SIM_MISS", config$missing_rates[1],
                                         "SIM_MISS")
rep_start <- as.integer(Sys.getenv("SIM_REP_START", "1"))
rep_end <- as.integer(Sys.getenv("SIM_REP_END", as.character(config$B_sim)))
reps <- rep_start:rep_end
cores <- v2_get_job_cores()

message("V2 continuous simulation: methods=", paste(methods, collapse = ","),
        " reps=", length(reps), " sigma=", sigma_u2, " miss=", miss_rate,
        " cores=", cores)
config$sigma_u2_grid <- sigma_u2
config$missing_rates <- miss_rate
paths <- parallel::mclapply(reps, function(replicate_id) {
  v2_run_continuous_one(replicate_id, sigma_u2, miss_rate, methods, config)
}, mc.cores = cores, mc.preschedule = FALSE, mc.silent = TRUE)
dir.create(file.path(config$results_root, config$run_label), recursive = TRUE, showWarnings = FALSE)
saveRDS(paths, file.path(config$results_root, config$run_label,
                         paste0("continuous_paths_",
                                v2_number_tag("sigma_u2", sigma_u2), "_",
                                v2_rate_tag(miss_rate), "_",
                                v2_clean_tag(paste(methods, collapse = "_")), ".rds")))
message("Completed continuous simulation.")
