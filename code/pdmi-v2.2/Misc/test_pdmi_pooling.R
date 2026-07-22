## Author: Jiaqi Bi

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "Misc/test_pdmi_pooling.R"
script_file <- gsub("~+~", " ", script_file, fixed = TRUE)
code_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = FALSE)
source(file.path(code_root, "Shared", "source_all.R"))

omega_names <- v22_omega_names()
make_fit <- function(x, scale) {
  V <- diag(scale, length(omega_names))
  dimnames(V) <- list(omega_names, omega_names)
  list(convergence = TRUE, omega = setNames(x, omega_names), vcov_omega = V)
}

fits <- list(
  make_fit(c(1, 2, 3, 4, 5), 0.1),
  make_fit(c(2, 4, 6, 8, 10), 0.2),
  make_fit(c(3, 6, 9, 12, 15), 0.3)
)
pool <- v22_pool_rubin_omega(fits, M_imp = 3L, require_all = TRUE)
mat <- do.call(rbind, lapply(fits, function(f) f$omega[omega_names]))
W_expected <- diag(mean(c(0.1, 0.2, 0.3)), length(omega_names))
dimnames(W_expected) <- list(omega_names, omega_names)
B_expected <- stats::cov(mat)
T_expected <- W_expected + (1 + 1 / 3) * B_expected
stopifnot(max(abs(pool$omega - colMeans(mat))) < 1e-12)
stopifnot(max(abs(pool$W - W_expected)) < 1e-12)
stopifnot(max(abs(pool$B - B_expected)) < 1e-12)
stopifnot(max(abs(pool$vcov_omega - T_expected)) < 1e-12)

bad_fits <- fits
bad_fits[[2]]$convergence <- FALSE
err <- try(v22_pool_rubin_omega(bad_fits, M_imp = 3L, require_all = TRUE), silent = TRUE)
stopifnot(inherits(err, "try-error"))

cat("V2.2 PDMI Rubin pooling tests passed\n")
