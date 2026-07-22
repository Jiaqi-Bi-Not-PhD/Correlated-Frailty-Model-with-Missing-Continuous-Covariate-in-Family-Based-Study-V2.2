## Author: Jiaqi Bi

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "Misc/test_exact_slice_pdmi.R"
script_file <- gsub("~+~", " ", script_file, fixed = TRUE)
code_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = FALSE)
source(file.path(code_root, "Shared", "source_all.R"))

cfg <- v22_default_config()
cfg$theta_slice_sweeps <- 1L
cfg$theta_slice_m <- 10L
cfg$ess_max_shrink <- 50L
cfg$gh_order <- 10L
omega <- v22_actual_omega(0.2, cfg)

dat <- data.frame(
  t0 = c(0, 0, 0, 0),
  time = c(45, 52, 46, 58),
  status = c(1, 0, 1, 0),
  mgene = c(1, 0, 1, 0),
  newx = c(0.1, -0.2, 0.3, 0.0),
  famID = c(1, 1, 2, 2),
  proband = c(1, 0, 1, 0),
  currentage = c(45, 52, 46, 58),
  indID = paste0("id", 1:4)
)
K <- diag(4)
rownames(K) <- colnames(K) <- dat$indID
U <- rep(0, nrow(dat))

stopifnot(is.finite(v22_log_prior_omega_exact(omega, cfg)))
bad <- omega
bad["sigma_u2"] <- -0.1
stopifnot(!is.finite(v22_log_prior_omega_exact(bad, cfg)))

lp <- v22_log_augmented_disease_posterior(dat, K, omega, U, cfg)
stopifnot(is.finite(lp))
lp_marginal <- v22_log_marginal_disease_posterior(dat, K, omega, cfg)
stopifnot(is.finite(lp_marginal))

set.seed(42)
draw <- v22_draw_theta_exact_slice(dat, K, omega, U, cfg)
stopifnot(all(is.finite(draw[v22_omega_names()])))
stopifnot(draw["sigma_u2"] > 0)
stopifnot(isTRUE(attr(draw, "exact_slice_diagnostics")$sigma_u2_sampled_directly))
stopifnot(isFALSE(attr(draw, "exact_slice_diagnostics")$auxiliary_u_used_in_theta_target))
stopifnot(identical(attr(draw, "exact_slice_diagnostics")$target,
                    "frailty_marginal_selected_likelihood_laplace_logscale_prior_no_jacobian"))
stopifnot(isFALSE(attr(draw, "exact_slice_diagnostics")$theta_prior_jacobian_log_rho_log_lambda))
stopifnot(isTRUE(attr(draw, "exact_slice_diagnostics")$no_fit_b_or_varHtotal_draw))

src_cont <- paste(deparse(body(v22_draw_continuous_pdmi)), collapse = "\n")
src_run_cont <- paste(deparse(body(v22_run_continuous_pdmi)), collapse = "\n")
stopifnot(!grepl("draw_omega_posterior", src_cont, fixed = TRUE))
stopifnot(!grepl("fit_frailtypack\\(", src_cont))
stopifnot(!grepl("fit_mean_completed_initial", src_run_cont, fixed = TRUE))

cat("V2.2 continuous exact-slice PDMI unit checks passed\n")
