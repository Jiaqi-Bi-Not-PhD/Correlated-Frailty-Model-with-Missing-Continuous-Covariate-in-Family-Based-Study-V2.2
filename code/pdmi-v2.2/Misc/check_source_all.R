## Author: Jiaqi Bi

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "Misc/check_source_all.R"
script_file <- gsub("~+~", " ", script_file, fixed = TRUE)
code_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = FALSE)
source(file.path(code_root, "Shared", "source_all.R"))

cfg <- v22_default_config()
stopifnot(identical(as.integer(cfg$M_imp_pdmi), 20L))
stopifnot(all(cfg$penetrance_prs == c(-0.5, 0, 0.5)))
stopifnot(identical(v22_method_tag("C-O-PDMI", "C-O"), "c_o_pdmi"))
stopifnot(identical(v22_method_tag("C-R-PDMI", "C-R"), "c_r_pdmi"))

required <- c(
  "v22_run_continuous_pdmi",
  "v22_pool_rubin_omega",
  "v22_penetrance_from_imputed_fits"
)
missing <- required[!vapply(required, exists, logical(1), mode = "function")]
if (length(missing)) stop("Missing required V2.2 functions: ", paste(missing, collapse = ", "))

cat("V2.2 continuous-covariate source check passed\n")
