## Author: Jiaqi Bi

## Lightweight structural check: source all modules without running simulations.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "Misc/check_source_all.R"
script_file <- gsub("~+~", " ", script_file, fixed = TRUE)
code_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = FALSE)
source(file.path(code_root, "Shared", "source_all.R"))

required_functions <- c(
  "v2_generate_complete_data",
  "v2_apply_continuous_missingness",
  "v2_run_cca",
  "v2_run_smcfcs_comparator",
  "v2_fit_frailtypack",
  "v2_penetrance_grid"
)
missing <- required_functions[!vapply(required_functions, exists, logical(1), mode = "function")]
if (length(missing)) stop("Missing required function(s): ", paste(missing, collapse = ", "))
message(
  "Continuous comparator source check passed: ",
  length(required_functions),
  " required functions found."
)
