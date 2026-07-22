## Author: Jiaqi Bi

## ============================================================
## frailtypack complete-covariate disease analysis.
## The final disease likelihood conditions on completed covariates and
## uses frailtypack's ascertainment-corrected correlated frailty fit.
## ============================================================

v2_fit_failure <- function(reason, method = NA_character_) {
  omega <- setNames(rep(NA_real_, 5), v2_omega_names())
  V <- matrix(NA_real_, 5, 5, dimnames = list(v2_omega_names(), v2_omega_names()))
  list(
    omega = omega,
    vcov_omega = V,
    convergence = FALSE,
    convergence_code = NA_integer_,
    failure_reason = reason,
    method = method,
    fit = NULL
  )
}

v2_extract_frailtypack_omega <- function(fit) {
  c(
    log.rho = log(as.numeric(fit$shape.weib[1])),
    log.lambda = log(as.numeric(fit$scale.weib[1])),
    beta_b = unname(fit$coef["mgene"]),
    beta_c = unname(fit$coef["newx"]),
    sigma_u2 = as.numeric(fit$sigma2)
  )
}

v2_match_frailtypack_coef <- function(coef_names, candidates, label) {
  hit <- match(candidates, coef_names)
  hit <- hit[!is.na(hit)]
  if (!length(hit)) {
    stop("Could not match ", label, " in fit$coef names: ",
         paste(coef_names, collapse = ", "))
  }
  hit[[1]]
}

v2_transform_frailtypack_covariance <- function(fit, coef_map = NULL) {
  raw <- as.numeric(fit$b)
  Vraw <- tryCatch(as.matrix(fit$varHtotal), error = function(e) NULL)
  if (!is.numeric(raw) || is.null(Vraw) || !is.matrix(Vraw)) {
    stop("frailtypack fit is missing fit$b or fit$varHtotal.")
  }
  if (length(raw) != nrow(Vraw) || nrow(Vraw) != ncol(Vraw)) {
    stop("fit$b and fit$varHtotal dimensions do not agree.")
  }
  if (any(!is.finite(raw)) || any(!is.finite(Vraw))) {
    stop("fit$b or fit$varHtotal contains non-finite values.")
  }

  p <- length(fit$coef)
  np <- length(raw)
  if (p < 2L || np <= p) stop("Unexpected frailtypack parameter layout.")

  ## frailtypack 2.13 optimizes positive Weibull and frailty-variance
  ## parameters on a signed square-root scale:
  ##   rho = b_rho^2, lambda = b_lambda^2, sigma_u2 = b_sigma^2.
  ## Use the signed b values so cross-covariances are transformed correctly.
  idx_shape_raw <- 1L
  idx_scale_raw <- 2L
  idx_sigma_raw <- np - p
  idx_coef <- seq.int(np - p + 1L, np)

  if (raw[idx_shape_raw] == 0 || raw[idx_scale_raw] == 0) {
    stop("Cannot transform Weibull covariance because a raw square-root parameter is zero.")
  }

  coef_names <- names(fit$coef)
  if (is.null(coef_names) || length(coef_names) != p || any(!nzchar(coef_names))) {
    stop("fit$coef must have names for covariance transformation.")
  }
  beta_b_candidates <- c(coef_map$beta_b %||% character(0),
                         "mgene", "beta_mgene", "beta_b", "majorgene")
  beta_c_candidates <- c(coef_map$beta_c %||% character(0),
                         "newx", "beta_PRS", "beta_c", "PRS")
  idx_beta_b <- idx_coef[v2_match_frailtypack_coef(coef_names, beta_b_candidates, "beta_b")]
  idx_beta_c <- idx_coef[v2_match_frailtypack_coef(coef_names, beta_c_candidates, "beta_c")]

  J <- matrix(0, nrow = length(v2_omega_names()), ncol = np)
  rownames(J) <- v2_omega_names()
  colnames(J) <- paste0("b", seq_len(np))
  colnames(J)[c(idx_shape_raw, idx_scale_raw, idx_sigma_raw,
                idx_beta_b, idx_beta_c)] <- c("sqrt.rho", "sqrt.lambda",
                                              "sqrt.sigma_u2", "beta_b", "beta_c")

  J["log.rho", idx_shape_raw] <- 2 / raw[idx_shape_raw]
  J["log.lambda", idx_scale_raw] <- 2 / raw[idx_scale_raw]
  J["beta_b", idx_beta_b] <- 1
  J["beta_c", idx_beta_c] <- 1
  J["sigma_u2", idx_sigma_raw] <- 2 * raw[idx_sigma_raw]

  V_reported <- J %*% Vraw %*% t(J)
  V_reported <- 0.5 * (V_reported + t(V_reported))
  dimnames(V_reported) <- list(v2_omega_names(), v2_omega_names())

  list(
    theta = v2_extract_frailtypack_omega(fit)[v2_omega_names()],
    vcov = V_reported[v2_omega_names(), v2_omega_names(), drop = FALSE],
    jacobian = J
  )
}

v2_extract_frailtypack_vcov <- function(fit) {
  nm <- v2_omega_names()
  transformed <- tryCatch(v2_transform_frailtypack_covariance(fit),
                          error = function(e) NULL)
  if (is.null(transformed) || any(!is.finite(transformed$vcov))) {
    return(matrix(NA_real_, 5, 5, dimnames = list(nm, nm)))
  }
  transformed$vcov[nm, nm, drop = FALSE]
}

v2_validate_analysis_dat <- function(dat, K) {
  if (!all(c("t0", "time", "status", "mgene", "newx", "famID",
             "proband", "currentage", "indID") %in% names(dat))) {
    stop("Analysis data are missing required frailtypack columns.")
  }
  if (anyNA(dat[, c("time", "status", "mgene", "newx", "famID", "proband", "currentage")])) {
    stop("Completed analysis data contain NA values.")
  }
  if (!v2_check_popplus_support(dat)) stop("Analysis data violate pop+ support.")
  K <- v2_align_K(K, dat)
  blocks <- v2_family_blocks(dat)
  for (idx in blocks) {
    Ki <- K[idx, idx, drop = FALSE]
    if (!v2_is_psd(Ki, tol = 1e-6)) stop("A family K_i is not positive semidefinite.")
    tryCatch(v2_safe_chol(Ki), error = function(e) stop("A family K_i is not positive definite."))
  }
  invisible(TRUE)
}

v2_fit_frailtypack <- function(dat, K, config = v2_default_config(), init_omega = NULL,
                               method = "frailtypack") {
  v2_require_packages(include_frailtypack = TRUE)
  dat <- as.data.frame(dat)
  if (!"t0" %in% names(dat)) dat$t0 <- 0
  K <- v2_align_K(K, dat)
  check <- tryCatch(v2_validate_analysis_dat(dat, K), error = function(e) e)
  if (inherits(check, "error")) return(v2_fit_failure(conditionMessage(check), method))

  init_omega <- init_omega %||% v2_actual_omega(config$sigma_u2_grid[1], config)
  init_b <- c(beta_b = unname(init_omega["beta_b"]),
              beta_c = unname(init_omega["beta_c"]))

  cluster <- survival::cluster
  Surv <- survival::Surv
  fit <- tryCatch({
    null <- if (.Platform$OS.type == "windows") "NUL" else "/dev/null"
    con <- file(null, open = "wt")
    sink(con)
    sink(con, type = "message")
    on.exit({
      sink(type = "message")
      sink()
      close(con)
    }, add = TRUE)
    suppressWarnings(
      frailtypack::frailtyPenal(
        Surv(t0, time, status) ~ mgene + newx + cluster(famID),
        data = dat,
        hazard = "Weibull",
        RandDist = "LogN",
        print.times = FALSE,
        init.B = init_b,
        covMatrix1 = as.matrix(K),
        recurrentAG = TRUE,
        maxit = config$frailtypack_maxit,
        proband = dat$proband,
        currentage = dat$currentage,
        agemin = config$agemin
      )
    )
  }, error = function(e) e)
  if (inherits(fit, "error")) return(v2_fit_failure(conditionMessage(fit), method))

  omega <- tryCatch(v2_extract_frailtypack_omega(fit), error = function(e) NULL)
  if (is.null(omega) || any(!is.finite(omega)) || omega["sigma_u2"] <= 0) {
    return(v2_fit_failure("Non-finite frailtypack estimate.", method))
  }
  V <- v2_extract_frailtypack_vcov(fit)
  if (any(!is.finite(V))) {
    return(v2_fit_failure("Non-finite frailtypack covariance on reported scale.", method))
  }
  if (!v2_is_psd(V)) V <- v2_near_psd(V)
  list(
    omega = omega[v2_omega_names()],
    vcov_omega = V[v2_omega_names(), v2_omega_names(), drop = FALSE],
    convergence = TRUE,
    convergence_code = fit$istop %||% NA_integer_,
    failure_reason = NA_character_,
    method = method,
    fit = fit
  )
}

v2_fit_mean_completed_initial <- function(dat, K, config) {
  tmp <- dat
  if (anyNA(tmp$newx)) {
    obs_mean <- mean(tmp$newx, na.rm = TRUE)
    if (!is.finite(obs_mean)) obs_mean <- 0
    tmp$newx[is.na(tmp$newx)] <- obs_mean
  }
  v2_fit_frailtypack(tmp, K, config, method = "initial_mean_completed")
}

v2_fit_list_to_psi <- function(disease_fit, nuisance_fit = NULL) {
  if (!isTRUE(disease_fit$convergence)) return(NULL)
  omega <- disease_fit$omega[v2_omega_names()]
  Vd <- disease_fit$vcov_omega[v2_omega_names(), v2_omega_names(), drop = FALSE]
  if (is.null(nuisance_fit) || length(nuisance_fit$eta) == 0L) {
    return(list(psi = omega, W = Vd))
  }
  eta <- nuisance_fit$eta
  Veta <- nuisance_fit$vcov
  nms <- c(names(omega), names(eta))
  W <- matrix(0, length(nms), length(nms), dimnames = list(nms, nms))
  W[names(omega), names(omega)] <- Vd
  W[names(eta), names(eta)] <- Veta
  list(psi = c(omega, eta), W = W)
}

v2_pool_mlmi_fullspace <- function(psi_fits, disease_names = v2_omega_names(), M_imp,
                                   eps = 1e-6) {
  ok <- vapply(psi_fits, function(x) is.list(x) && !is.null(x$psi) &&
                 all(is.finite(x$psi)) && !is.null(x$W), logical(1))
  psi_fits <- psi_fits[ok]
  if (!length(psi_fits)) stop("No successful completed-data fits to pool.")
  common_names <- names(psi_fits[[1]]$psi)
  par_mat <- do.call(rbind, lapply(psi_fits, function(x) x$psi[common_names]))
  M <- nrow(par_mat)
  psi_bar <- colMeans(par_mat)
  W <- Reduce("+", lapply(psi_fits, function(x) {
    V <- x$W[common_names, common_names, drop = FALSE]
    V[!is.finite(V)] <- 0
    V
  })) / M
  B <- if (M > 1L) stats::cov(par_mat) else matrix(0, length(common_names), length(common_names))
  dimnames(B) <- list(common_names, common_names)
  W <- v2_near_psd(W, eps)
  B <- v2_near_psd(B, eps = 0)

  Winv_half <- v2_matrix_sqrt_sym(W, inverse = TRUE, eps = eps)
  Whalf <- v2_matrix_sqrt_sym(W, inverse = FALSE, eps = eps)
  A <- 0.5 * (Winv_half %*% B %*% Winv_half + t(Winv_half %*% B %*% Winv_half))
  ee <- eigen(A, symmetric = TRUE)
  lam <- pmin(pmax(ee$values, 0), 1 - eps)
  A_trunc <- ee$vectors %*% diag(lam, length(lam)) %*% t(ee$vectors)
  V_ML <- Whalf %*% solve(diag(ncol(W)) - A_trunc) %*% Whalf
  T_MLMI <- V_ML + B / M_imp
  dimnames(V_ML) <- dimnames(T_MLMI) <- list(common_names, common_names)
  list(
    psi_bar = psi_bar,
    omega = psi_bar[disease_names],
    V_ML_full = V_ML,
    T_MLMI_full = T_MLMI,
    vcov_omega = T_MLMI[disease_names, disease_names, drop = FALSE],
    W = W,
    B = B,
    gamma_mis_eigen = lam,
    M_success = M
  )
}

v2_pool_rubin_omega <- function(fits, M_imp) {
  ok <- vapply(fits, function(x) isTRUE(x$convergence) && all(is.finite(x$omega)), logical(1))
  fits <- fits[ok]
  if (!length(fits)) stop("No successful completed-data fits for Rubin pooling.")
  omega_mat <- do.call(rbind, lapply(fits, function(x) x$omega[v2_omega_names()]))
  M <- nrow(omega_mat)
  qbar <- colMeans(omega_mat)
  W <- Reduce("+", lapply(fits, function(x) x$vcov_omega[v2_omega_names(), v2_omega_names(), drop = FALSE])) / M
  B <- if (M > 1L) stats::cov(omega_mat) else matrix(0, 5, 5)
  dimnames(B) <- list(v2_omega_names(), v2_omega_names())
  T <- W + (1 + 1 / M_imp) * B
  list(omega = qbar, vcov_omega = v2_near_psd(T), W = W, B = B, M_success = M)
}

v2_penetrance_from_fit <- function(omega, vcov_omega, config = v2_default_config()) {
  grid <- expand.grid(
    age = config$penetrance_ages,
    prs = config$penetrance_prs,
    gene = config$penetrance_gene,
    KEEP.OUT.ATTRS = FALSE
  )
  grid$estimate <- mapply(v2_penetrance, grid$age, grid$prs, grid$gene,
                          MoreArgs = list(omega = omega, k0 = config$penetrance_k0,
                                          agemin = config$agemin,
                                          gh_order = config$gh_order))
  grid$se <- mapply(function(age, prs, gene) {
    g <- v2_penetrance_gradient_omega(age, prs, gene, omega,
                                      k0 = config$penetrance_k0,
                                      agemin = config$agemin,
                                      gh_order = config$gh_order)
    sqrt(max(as.numeric(t(g) %*% vcov_omega[v2_omega_names(), v2_omega_names()] %*% g), 0))
  }, grid$age, grid$prs, grid$gene)
  grid
}

v2_penetrance_from_imputed_fits <- function(fits, pooled_vcov_omega, M_imp,
                                            config = v2_default_config()) {
  ok <- vapply(fits, function(x) isTRUE(x$convergence) && all(is.finite(x$omega)), logical(1))
  fits <- fits[ok]
  if (!length(fits)) stop("No successful fits for penetrance.")
  grid0 <- expand.grid(
    age = config$penetrance_ages,
    prs = config$penetrance_prs,
    gene = config$penetrance_gene,
    KEEP.OUT.ATTRS = FALSE
  )
  rows <- lapply(seq_len(nrow(grid0)), function(r) {
    qmat <- do.call(rbind, lapply(fits, function(f) {
      v2_penetrance(grid0$age[r], grid0$prs[r], grid0$gene[r], f$omega,
                    k0 = config$penetrance_k0, agemin = config$agemin,
                    gh_order = config$gh_order)
    }))
    est <- mean(qmat)
    Bq <- if (length(qmat) > 1L) stats::var(as.numeric(qmat)) else 0
    grad <- Reduce("+", lapply(fits, function(f) {
      v2_penetrance_gradient_omega(grid0$age[r], grid0$prs[r], grid0$gene[r],
                                   f$omega, k0 = config$penetrance_k0,
                                   agemin = config$agemin,
                                   gh_order = config$gh_order)
    })) / length(fits)
    var <- as.numeric(t(grad) %*% pooled_vcov_omega[v2_omega_names(), v2_omega_names()] %*% grad) +
      Bq / M_imp
    data.frame(grid0[r, , drop = FALSE], estimate = est, se = sqrt(max(var, 0)))
  })
  do.call(rbind, rows)
}

v2_penetrance_from_mlmi_streams <- function(point_fits, V_ML_omega, M_point,
                                           variance_fits = NULL,
                                           config = v2_default_config()) {
  ok_point <- vapply(point_fits, function(x) isTRUE(x$convergence) && all(is.finite(x$omega)), logical(1))
  point_fits <- point_fits[ok_point]
  if (!length(point_fits)) stop("No successful point-estimation fits for penetrance.")
  variance_fits <- variance_fits %||% point_fits
  ok_var <- vapply(variance_fits, function(x) isTRUE(x$convergence) && all(is.finite(x$omega)), logical(1))
  variance_fits <- variance_fits[ok_var]
  if (!length(variance_fits)) stop("No successful variance-stream fits for penetrance.")
  grid0 <- expand.grid(
    age = config$penetrance_ages,
    prs = config$penetrance_prs,
    gene = config$penetrance_gene,
    KEEP.OUT.ATTRS = FALSE
  )
  rows <- lapply(seq_len(nrow(grid0)), function(r) {
    q_point <- do.call(rbind, lapply(point_fits, function(f) {
      v2_penetrance(grid0$age[r], grid0$prs[r], grid0$gene[r], f$omega,
                    k0 = config$penetrance_k0, agemin = config$agemin,
                    gh_order = config$gh_order)
    }))
    q_var <- do.call(rbind, lapply(variance_fits, function(f) {
      v2_penetrance(grid0$age[r], grid0$prs[r], grid0$gene[r], f$omega,
                    k0 = config$penetrance_k0, agemin = config$agemin,
                    gh_order = config$gh_order)
    }))
    est <- mean(q_point)
    Bq <- if (length(q_var) > 1L) stats::var(as.numeric(q_var)) else 0
    grad <- Reduce("+", lapply(point_fits, function(f) {
      v2_penetrance_gradient_omega(grid0$age[r], grid0$prs[r], grid0$gene[r],
                                   f$omega, k0 = config$penetrance_k0,
                                   agemin = config$agemin,
                                   gh_order = config$gh_order)
    })) / length(point_fits)
    var <- as.numeric(t(grad) %*% V_ML_omega[v2_omega_names(), v2_omega_names()] %*% grad) +
      Bq / M_point
    data.frame(grid0[r, , drop = FALSE], estimate = est, se = sqrt(max(var, 0)))
  })
  do.call(rbind, rows)
}
