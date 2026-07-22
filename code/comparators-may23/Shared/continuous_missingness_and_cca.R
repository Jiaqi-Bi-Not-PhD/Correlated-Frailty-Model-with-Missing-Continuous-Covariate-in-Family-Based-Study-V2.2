## Author: Jiaqi Bi

## ============================================================
## Continuous-PRS MAR masks and CCA row-aligned subsetting.
## Missingness is post-generation and fixed-size weighted.
## ============================================================

v2_family_size_z <- function(dat) {
  n_i <- ave(dat$indID, dat$famID, FUN = length)
  s <- stats::sd(unique(n_i))
  if (!is.finite(s) || s <= 0) return(rep(0, nrow(dat)))
  (n_i - mean(unique(n_i))) / s
}

v2_weighted_fixed_mask <- function(weights, n_miss) {
  n_miss <- as.integer(n_miss)
  miss <- rep(FALSE, length(weights))
  if (n_miss <= 0L) return(miss)
  if (n_miss >= length(weights)) {
    miss[] <- TRUE
    return(miss)
  }
  weights <- as.numeric(weights)
  weights[!is.finite(weights) | weights < 0] <- 0
  if (!any(weights > 0)) weights <- rep(1, length(weights))
  pick <- sample(seq_along(weights), size = n_miss, replace = FALSE, prob = weights)
  miss[pick] <- TRUE
  miss
}

v2_apply_continuous_missingness <- function(dat, target_missing_rate, seed, config = v2_default_config()) {
  set.seed(seed)
  out <- dat
  fz <- v2_family_size_z(out)
  score <- 0.35 * as.numeric(out$proband) +
    0.30 * as.numeric(out$status) +
    0.25 * as.numeric(out$mgene) -
    0.25 * (as.numeric(out$time) / 100) +
    0.15 * fz
  eligible <- rep(TRUE, nrow(out))
  n_miss <- round(target_missing_rate * sum(eligible))
  miss_eligible <- v2_weighted_fixed_mask(exp(pmin(score[eligible], 700)), n_miss)
  miss <- rep(FALSE, nrow(out))
  miss[eligible] <- miss_eligible
  out$newx[miss] <- NA_real_
  attr(out, "missing_mask") <- miss
  attr(out, "missing_score") <- score
  attr(out, "missing_type") <- "continuous"
  out
}
v2_missing_diagnostics <- function(dat_with_missing, target_missing_rate) {
  mask <- attr(dat_with_missing, "missing_mask")
  if (is.null(mask)) mask <- is.na(dat_with_missing$newx)
  incomplete_families <- length(unique(dat_with_missing$famID[mask]))
  data.frame(
    target_missing_rate = target_missing_rate,
    eligible_n = nrow(dat_with_missing),
    missing_n = sum(mask),
    realized_missing_rate = sum(mask) / nrow(dat_with_missing),
    missing_proband_prs_count = sum(
      mask & as.numeric(dat_with_missing$proband) == 1
    ),
    incomplete_family_count = incomplete_families,
    stringsAsFactors = FALSE
  )
}

v2_make_continuous_cca_dataset <- function(dat, K) {
  dat <- as.data.frame(dat)
  K <- v2_align_K(K, dat)
  keep <- rep(TRUE, nrow(dat))
  removed_families <- character(0)
  blocks <- v2_family_blocks(dat)
  for (idx in blocks) {
    proband_index <- idx[as.numeric(dat$proband[idx]) == 1]
    if (length(proband_index) != 1L) {
      keep[idx] <- FALSE
      removed_families <- c(
        removed_families,
        as.character(dat$famID[idx[1]])
      )
    } else if (is.na(dat$newx[proband_index])) {
      keep[idx] <- FALSE
      removed_families <- c(
        removed_families,
        as.character(dat$famID[idx[1]])
      )
    } else {
      keep[idx[is.na(dat$newx[idx])]] <- FALSE
    }
  }
  dat_cc <- dat[keep, , drop = FALSE]
  K_cc <- K[keep, keep, drop = FALSE]
  rownames(K_cc) <- colnames(K_cc) <- as.character(dat_cc$indID)
  diagnostics <- data.frame(
    n_original = nrow(dat),
    n_retained = nrow(dat_cc),
    n_removed = nrow(dat) - nrow(dat_cc),
    families_original = length(unique(dat$famID)),
    families_retained = length(unique(dat_cc$famID)),
    families_removed_whole = length(unique(removed_families)),
    stringsAsFactors = FALSE
  )
  list(
    dat = dat_cc,
    K = v2_align_K(K_cc, dat_cc),
    diagnostics = diagnostics
  )
}

v2_run_cca <- function(dat_with_missing, K, config = v2_default_config(),
                       init_omega = NULL) {
  cc <- tryCatch(
    v2_make_continuous_cca_dataset(dat_with_missing, K),
    error = function(e) e
  )
  if (inherits(cc, "error")) {
    fit <- v2_fit_failure(conditionMessage(cc), "CCA")
    return(list(fit = fit, diagnostics = data.frame()))
  }
  fit <- v2_fit_frailtypack(
    cc$dat,
    cc$K,
    config,
    init_omega,
    method = "CCA"
  )
  list(fit = fit, diagnostics = cc$diagnostics)
}
