## Author: Jiaqi Bi

## Source the continuous-covariate comparator modules in dependency order.

v2_this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
v2_this_file <- gsub("~+~", " ", v2_this_file, fixed = TRUE)
v2_this_file <- tryCatch(normalizePath(v2_this_file), error = function(e) NA_character_)
if (is.na(v2_this_file) || !nzchar(v2_this_file)) {
  wd <- normalizePath(getwd(), mustWork = FALSE)
  v2_code_root <- if (basename(wd) == "Shared") dirname(wd) else wd
} else {
  v2_code_root <- normalizePath(file.path(dirname(v2_this_file), ".."))
}

source(file.path(v2_code_root, "Shared", "config.R"))
source(file.path(v2_code_root, "Shared", "packages_and_sources.R"))
source(file.path(v2_code_root, "Shared", "numerics.R"))
source(file.path(v2_code_root, "Shared", "data_generation.R"))
source(file.path(v2_code_root, "Shared", "frailtypack_analysis.R"))
source(file.path(v2_code_root, "Shared", "continuous_missingness_and_cca.R"))
source(file.path(v2_code_root, "Shared", "smcfcs_continuous_comparator.R"))
source(file.path(v2_code_root, "Shared", "results_and_summaries.R"))

v2_set_thread_env()
