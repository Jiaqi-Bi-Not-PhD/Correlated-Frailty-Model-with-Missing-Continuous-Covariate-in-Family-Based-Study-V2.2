## Author: Jiaqi Bi

## ============================================================
## MI-SMCFCS comparator for a missing continuous PRS.
## This is intentionally not frailtypack-congenial: it does not use
## M_i(Y|X) or -log alpha_i^F for a missing proband PRS.
## Final completed-data fits still use frailtypack.
## ============================================================

v2_extract_smcfcs_impdatasets <- function(imp) {
  if (is.null(imp)) return(NULL)
  if (!is.null(imp$impDatasets)) return(imp$impDatasets)
  if (!is.null(imp$imp.dataset)) return(imp$imp.dataset)
  if (is.list(imp) && all(vapply(imp, is.data.frame, logical(1)))) return(imp)
  NULL
}

v2_run_smcfcs_imputation <- function(dat, M = 10L, seed = NULL,
                                     config = v2_default_config()) {
  if (!requireNamespace("smcfcs", quietly = TRUE)) {
    stop("Package smcfcs is not installed; MI-SMCFCS comparator cannot run.")
  }
  if (!is.null(seed)) set.seed(seed)
  smdat <- dat[
    ,
    c("time", "status", "mgene", "newx", "proband", "fsize"),
    drop = FALSE
  ]
  smdat$status <- as.integer(smdat$status)
  smdat$proband <- as.integer(smdat$proband)
  smdat$fsize <- as.numeric(smdat$fsize)

  method <- rep("", ncol(smdat))
  names(method) <- names(smdat)
  method["newx"] <- "norm"

  pred <- matrix(
    0,
    ncol(smdat),
    ncol(smdat),
    dimnames = list(names(smdat), names(smdat))
  )
  pred["newx", c("time", "status", "mgene", "proband", "fsize")] <- 1
  diag(pred) <- 0

  imp <- tryCatch(
    smcfcs::smcfcs(
      originaldata = smdat,
      smtype = "weibull",
      smformula = "Surv(time, status) ~ mgene + newx",
      method = method,
      predictorMatrix = pred,
      m = M,
      numit = 20L,
      rjlimit = 5000,
      noisy = FALSE
    ),
    error = function(e) e
  )
  if (inherits(imp, "error")) stop(conditionMessage(imp))

  imp_list <- v2_extract_smcfcs_impdatasets(imp)
  if (is.null(imp_list) || !length(imp_list)) {
    stop("smcfcs returned no completed datasets.")
  }
  completed <- lapply(imp_list, function(x) {
    di <- dat
    di$newx <- as.numeric(x$newx)
    di
  })
  list(completed = completed, smcfcs = imp)
}

v2_run_smcfcs_comparator <- function(dat_with_missing, K,
                                     config = v2_default_config(),
                                     seed = NULL) {
  imp <- tryCatch(
    v2_run_smcfcs_imputation(
      dat_with_missing,
      M = config$M_imp_smcfcs,
      seed = seed,
      config = config
    ),
    error = function(e) e
  )
  if (inherits(imp, "error")) {
    return(list(
      convergence = FALSE,
      failure_reason = conditionMessage(imp),
      imputations = NULL
    ))
  }

  init_fit <- v2_fit_mean_completed_initial(dat_with_missing, K, config)
  init_omega <- if (isTRUE(init_fit$convergence)) init_fit$omega else NULL
  fits <- lapply(
    imp$completed,
    v2_fit_frailtypack,
    K = K,
    config = config,
    init_omega = init_omega,
    method = "MI-SMCFCS"
  )
  pool <- tryCatch(
    v2_pool_rubin_omega(fits, M_imp = config$M_imp_smcfcs),
    error = function(e) e
  )
  if (inherits(pool, "error")) {
    return(list(
      convergence = FALSE,
      failure_reason = conditionMessage(pool),
      disease_fits = fits,
      imputations = imp
    ))
  }
  pen <- v2_penetrance_from_imputed_fits(
    fits,
    pool$vcov_omega,
    M_imp = config$M_imp_smcfcs,
    config = config
  )
  list(
    convergence = TRUE,
    failure_reason = NA_character_,
    pooled = pool,
    penetrance = pen,
    disease_fits = fits,
    imputations = imp,
    diagnostics = list(m_success = pool$M_success)
  )
}
