## Author: Jiaqi Bi

## ============================================================
## Continuous-PRS MAR masks for the V2.2 PDMI simulations.
## Missingness is post-generation and fixed-size weighted.
## ============================================================

v22_family_size_z <- function(dat) {
  n_i <- ave(dat$indID, dat$famID, FUN = length)
  s <- stats::sd(unique(n_i))
  if (!is.finite(s) || s <= 0) return(rep(0, nrow(dat)))
  (n_i - mean(unique(n_i))) / s
}

v22_weighted_fixed_mask <- function(weights, n_miss) {
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

v22_apply_continuous_missingness <- function(dat, target_missing_rate, seed, config = v22_default_config()) {
  set.seed(seed)
  out <- dat
  fz <- v22_family_size_z(out)
  score <- 0.35 * as.numeric(out$proband) +
    0.30 * as.numeric(out$status) +
    0.25 * as.numeric(out$mgene) -
    0.25 * (as.numeric(out$time) / 100) +
    0.15 * fz
  eligible <- rep(TRUE, nrow(out))
  n_miss <- round(target_missing_rate * sum(eligible))
  miss_eligible <- v22_weighted_fixed_mask(exp(pmin(score[eligible], 700)), n_miss)
  miss <- rep(FALSE, nrow(out))
  miss[eligible] <- miss_eligible
  out$newx[miss] <- NA_real_
  attr(out, "missing_mask") <- miss
  attr(out, "missing_score") <- score
  attr(out, "missing_type") <- "continuous"
  out
}
v22_missing_diagnostics <- function(dat_with_missing, target_missing_rate) {
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
