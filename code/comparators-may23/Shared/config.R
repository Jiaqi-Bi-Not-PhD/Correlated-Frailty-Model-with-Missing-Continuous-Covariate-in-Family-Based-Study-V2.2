## Author: Jiaqi Bi

## ============================================================
## V2 frailtypack-congenial simulation configuration
## Technical-note source: Simulation design, scenario grid, and
## reported scale omega = (log rho, log lambda, beta_b, beta_c, sigma_u^2).
## ============================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b

v2_omega_names <- function() {
  c("log.rho", "log.lambda", "beta_b", "beta_c", "sigma_u2")
}

v2_default_config <- function() {
  list(
    B_sim = 1000L,
    M_imp_smcfcs = 20L,
    n_families = 498L,
    sigma_u2_grid = c(0.5, 0.2),
    missing_rates = c(0.20, 0.50, 0.80),
    pm = 0.02,
    omega_base = c(log.rho = 0.804, log.lambda = 4.71,
                   beta_b = 2.2, beta_c = 1.0),
    design = "pop+",
    variation = "kinship",
    base.dist = "Weibull",
    frailty.dist = "lognormal",
    interaction = FALSE,
    probandage = c(45, 2),
    agemin = 0,
    agemax = 100,
    align_female_only = TRUE,
    analysis_female_only = TRUE,
    prs_sigma2 = 0.1,
    penetrance_ages = c(40, 50, 60, 70, 80),
    penetrance_prs = c(-0.5, 0, 0.5),
    penetrance_gene = c(0, 1),
    penetrance_k0 = 1,
    frailtypack_maxit = 35L,
    gh_order = 20L,
    skip_existing_results = TRUE,
    skip_existing_benchmarks = TRUE,
    run_label = Sys.getenv("SIM_RUN_LABEL", paste0(Sys.info()[["user"]], "_", format(Sys.time(), "%Y%m%d_%H%M%S"))),
    results_root = Sys.getenv("SIM_RESULTS_ROOT", file.path("Results", "raw"))
  )
}

v2_actual_omega <- function(sigma_u2, config = v2_default_config()) {
  c(config$omega_base, sigma_u2 = sigma_u2)[v2_omega_names()]
}

v2_actual_value_legacy_names <- function(sigma_u2, config = v2_default_config()) {
  c(log.shape = unname(config$omega_base["log.rho"]),
    log.scale = unname(config$omega_base["log.lambda"]),
    beta_mgene = unname(config$omega_base["beta_b"]),
    beta_PRS = unname(config$omega_base["beta_c"]),
    sigma2 = sigma_u2)
}

v2_theta_from_omega <- function(omega) {
  omega <- omega[v2_omega_names()]
  list(
    log_rho = unname(omega["log.rho"]),
    rho = exp(unname(omega["log.rho"])),
    log_lambda = unname(omega["log.lambda"]),
    lambda = exp(unname(omega["log.lambda"])),
    beta_b = unname(omega["beta_b"]),
    beta_c = unname(omega["beta_c"]),
    sigma_u2 = unname(omega["sigma_u2"])
  )
}

v2_log_kappa_from_omega <- function(omega) {
  th <- v2_theta_from_omega(omega)
  -th$rho * th$log_lambda
}

v2_method_tag <- function(method, prior_version = NA_character_) {
  method <- tolower(method)
  prior_version <- tolower(prior_version %||% "")
  if (identical(method, "full-data")) return("full_data")
  if (identical(method, "cca")) return("cca")
  if (identical(method, "mi-smcfcs")) return("mi_smcfcs")
  gsub("[^a-z0-9]+", "_", paste(method, prior_version))
}

v2_parse_numeric_env <- function(name, default) {
  value <- Sys.getenv(name, paste(default, collapse = ","))
  pieces <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  pieces <- pieces[nzchar(pieces)]
  out <- suppressWarnings(as.numeric(pieces))
  if (!length(out) || anyNA(out)) {
    stop("Environment variable ", name, " must contain numeric value(s), got: ", value)
  }
  out
}

v2_parse_single_numeric_env <- function(name, default, label = name) {
  values <- v2_parse_numeric_env(name, default)
  if (length(values) != 1L) {
    stop(label, " must contain exactly one value per job. Got ",
         name, "=", paste(values, collapse = ","),
         ". Submit separate jobs for repeated grid values.")
  }
  values[[1]]
}

v2_number_tag <- function(prefix, x) {
  paste0(prefix, "_", gsub("\\.", "p", as.character(x)))
}

v2_rate_tag <- function(x) {
  if (is.na(x)) return("nomiss")
  sprintf("miss%02d", round(100 * x))
}

v2_clean_tag <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  gsub("_+", "_", x)
}

v2_find_code_root <- function(starts = c(getwd())) {
  for (start in starts) {
    if (!nzchar(start)) next
    cur <- normalizePath(start, mustWork = FALSE)
    if (!dir.exists(cur)) cur <- dirname(cur)
    repeat {
      if (dir.exists(file.path(cur, "Shared")) &&
          dir.exists(file.path(cur, "Continuous Only"))) {
        return(normalizePath(cur))
      }
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
  }
  stop("Could not locate the V2 code root.")
}

v2_get_job_cores <- function(default = 1L) {
  if (identical(.Platform$OS.type, "windows")) return(1L)
  value <- Sys.getenv("SIM_CORES", unset = "")
  if (nzchar(value)) {
    n <- suppressWarnings(as.integer(value))
    if (is.na(n) || n < 1L) {
      stop("SIM_CORES must be a positive integer; got: ", value)
    }
    return(n)
  }
  n <- suppressWarnings(parallel::detectCores())
  if (!is.na(n) && n > 0L) return(n)
  as.integer(default)
}

v2_seed_streams <- function(replicate_id, missing_type = "none", target_missing_rate = NA_real_,
                            method = "full-data", seed_base = 700000L) {
  type_offset <- switch(
    missing_type,
    none = 0L,
    continuous = 10000L,
    stop("Unsupported missingness type in the continuous release: ", missing_type)
  )
  rate_offset <- if (is.na(target_missing_rate)) 0L else as.integer(round(1000 * target_missing_rate))
  method_offset <- sum(utf8ToInt(method)) %% 10000L
  list(
    complete_seed = seed_base + replicate_id,
    missing_mask_seed = seed_base + 100000L + type_offset + rate_offset + replicate_id,
    method_seed = seed_base + 200000L + type_offset + rate_offset + method_offset + replicate_id
  )
}

v2_result_metadata <- function(replicate_id, missing_type, sigma_u2, target_missing_rate,
                               method, prior_version, M_imp, seeds, config) {
  list(
    replicate_id = replicate_id,
    missing_type = missing_type,
    sigma_u2 = sigma_u2,
    target_missing_rate = target_missing_rate,
    method = method,
    prior_version = prior_version,
    M_imp = M_imp,
    complete_data_seed = seeds$complete_seed,
    missing_mask_seed = seeds$missing_mask_seed,
    method_seed = seeds$method_seed,
    run_label = config$run_label,
    run_owner = Sys.info()[["user"]],
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
}
