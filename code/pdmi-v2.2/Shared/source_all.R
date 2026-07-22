## Author: Jiaqi Bi

## Source the continuous-covariate V2.2 modules in dependency order.

v22_this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
v22_this_file <- gsub("~+~", " ", v22_this_file, fixed = TRUE)
v22_this_file <- tryCatch(normalizePath(v22_this_file), error = function(e) NA_character_)
if (is.na(v22_this_file) || !nzchar(v22_this_file)) {
  wd <- normalizePath(getwd(), mustWork = FALSE)
  v22_code_root <- if (basename(wd) == "Shared") dirname(wd) else wd
} else {
  v22_code_root <- normalizePath(file.path(dirname(v22_this_file), ".."))
}

source(file.path(v22_code_root, "Shared", "config.R"))
source(file.path(v22_code_root, "Shared", "packages_and_sources.R"))
source(file.path(v22_code_root, "Shared", "numerics.R"))
source(file.path(v22_code_root, "Shared", "pdmi_diagnostics.R"))
source(file.path(v22_code_root, "Shared", "data_generation.R"))
source(file.path(v22_code_root, "Shared", "frailtypack_analysis.R"))
source(file.path(v22_code_root, "Shared", "continuous_missingness.R"))
source(file.path(v22_code_root, "Shared", "pdmi_continuous.R"))
source(file.path(v22_code_root, "Shared", "pdmi_exact_slice_continuous.R"))
source(file.path(v22_code_root, "Shared", "results_and_summaries.R"))

v22_set_thread_env()
