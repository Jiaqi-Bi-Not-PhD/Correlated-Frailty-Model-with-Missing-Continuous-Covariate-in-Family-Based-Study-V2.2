## Author: Jiaqi Bi

## Unit checks for finite-value PDMI diagnostics.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "Misc/test_pdmi_diagnostics.R"
script_file <- gsub("~+~", " ", script_file, fixed = TRUE)
code_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = FALSE)
if (!dir.exists(file.path(code_root, "Shared"))) code_root <- normalizePath(getwd(), mustWork = FALSE)
source(file.path(code_root, "Shared", "source_all.R"))

stopifnot_true <- function(x, msg) {
  if (!isTRUE(x)) stop(msg, call. = FALSE)
}

dat <- data.frame(
  t0 = 0,
  time = c(50, 60, 70),
  status = c(1, 0, 1),
  mgene = c(1, 0, 1),
  newx = c(0.1, -0.2, 0.3),
  famID = c(1, 1, 1),
  proband = c(1, 0, 0),
  currentage = c(45, 50, 55),
  indID = c("a", "b", "c"),
  stringsAsFactors = FALSE
)
K <- diag(3)
rownames(K) <- colnames(K) <- dat$indID

d0 <- v22_analysis_data_diagnostics(dat, K, context = list(test = "clean"))
stopifnot_true(!d0$has_bad_input, "Clean complete data should not be flagged.")
stopifnot_true(nrow(dat) == d0$n_rows, "Diagnostic validation must not drop rows.")

d_na <- dat
d_na$newx[2] <- NA_real_
diag_na <- v22_analysis_data_diagnostics(d_na, K, context = list(test = "na"))
stopifnot_true(diag_na$has_bad_input, "NA newx should be flagged.")
stopifnot_true("newx" %in% diag_na$bad_columns, "NA newx should identify the bad column.")

d_nan <- dat
d_nan$newx[2] <- NaN
diag_nan <- v22_analysis_data_diagnostics(d_nan, K, context = list(test = "nan"))
stopifnot_true(diag_nan$has_bad_input, "NaN newx should be flagged.")

d_inf <- dat
d_inf$newx[2] <- Inf
diag_inf <- v22_analysis_data_diagnostics(d_inf, K, context = list(test = "inf"))
stopifnot_true(diag_inf$has_bad_input, "Inf newx should be flagged.")
stopifnot_true(isTRUE(diag_inf$has_inf), "Inf newx should set has_inf.")

d_ninf <- dat
d_ninf$newx[2] <- -Inf
diag_ninf <- v22_analysis_data_diagnostics(d_ninf, K, context = list(test = "neg_inf"))
stopifnot_true(diag_ninf$has_bad_input, "-Inf newx should be flagged.")

err <- tryCatch(
  v22_stop_pdmi_diagnostic("diagnostic failure", list(stage = list(m = 1L, iter = 2L))),
  error = function(e) e
)
stopifnot_true(inherits(err, "v22_pdmi_diagnostic_error"), "Diagnostic errors should preserve class.")
stopifnot_true(err$diagnostics$stage$m == 1L && err$diagnostics$stage$iter == 2L,
               "Diagnostic errors should preserve stage metadata.")

message("PDMI diagnostics tests passed.")
