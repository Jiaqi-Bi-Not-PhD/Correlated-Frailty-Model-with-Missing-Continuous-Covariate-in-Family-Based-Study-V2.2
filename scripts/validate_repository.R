## Author: Jiaqi Bi

## Structural, configuration, and dependency checks for the R-code release.
## Add --strict-packages to require the package versions used for validation.

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_file <- if (length(file_arg)) {
  sub("^--file=", "", file_arg[1])
} else {
  "scripts/validate_repository.R"
}
script_file <- normalizePath(script_file, mustWork = TRUE)
repository_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
strict_packages <- "--strict-packages" %in% commandArgs(trailingOnly = TRUE)

assert_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}

required_paths <- c(
  "scripts/install_dependencies.R",
  "scripts/run_simulation_grid.R",
  "code/pdmi-v2.2/Continuous Only/run_continuous_pdmi_simulation.R",
  "code/pdmi-v2.2/Shared/continuous_missingness.R",
  "code/pdmi-v2.2/Shared/pdmi_continuous.R",
  "code/pdmi-v2.2/Shared/pdmi_exact_slice_continuous.R",
  "code/pdmi-v2.2/Shared/source_all.R",
  "code/pdmi-v2.2/Misc/summarize_results.R",
  "code/comparators-may23/Continuous Only/run_continuous_simulation.R",
  "code/comparators-may23/No Missing/run_no_missing_benchmark.R",
  "code/comparators-may23/Shared/continuous_missingness_and_cca.R",
  "code/comparators-may23/Shared/smcfcs_continuous_comparator.R",
  "code/comparators-may23/Shared/source_all.R",
  "code/comparators-may23/Misc/summarize_results.R"
)
missing_paths <- required_paths[
  !file.exists(file.path(repository_root, required_paths))
]
assert_true(
  !length(missing_paths),
  paste("Missing required R-code paths:", paste(missing_paths, collapse = ", "))
)

pdmi_config_env <- new.env(parent = globalenv())
sys.source(
  file.path(repository_root, "code", "pdmi-v2.2", "Shared", "config.R"),
  envir = pdmi_config_env
)
pdmi_config <- pdmi_config_env$v22_default_config()
assert_true(
  identical(pdmi_config$B_sim, 1000L),
  "V2.2 PDMI must default to 1,000 attempted replicates per cell."
)
assert_true(
  identical(pdmi_config$M_imp_pdmi, 20L),
  "V2.2 PDMI must default to M=20."
)

comparator_config_env <- new.env(parent = globalenv())
sys.source(
  file.path(repository_root, "code", "comparators-may23", "Shared", "config.R"),
  envir = comparator_config_env
)
comparator_config <- comparator_config_env$v2_default_config()
assert_true(
  identical(comparator_config$B_sim, 1000L),
  "Comparator simulations must default to 1,000 attempted replicates per cell."
)
assert_true(
  identical(comparator_config$M_imp_smcfcs, 20L),
  "MI-SMCFCS must default to M=20."
)

packages <- c(
  "Matrix", "parallel", "survival", "kinship2", "MASS", "truncnorm",
  "frailtypack", "smcfcs"
)
available <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
versions <- vapply(packages, function(package) {
  if (requireNamespace(package, quietly = TRUE)) {
    as.character(utils::packageVersion(package))
  } else {
    NA_character_
  }
}, character(1))
package_table <- data.frame(
  package = packages,
  available = unname(available),
  version = unname(versions),
  stringsAsFactors = FALSE
)

tested_versions <- c(frailtypack = "2.13", smcfcs = "2.0.2")
tested_available <- available[names(tested_versions)]
tested_match <- tested_available &
  versions[names(tested_versions)] == tested_versions

if (strict_packages && any(!available)) {
  stop(
    "Missing required R packages: ",
    paste(packages[!available], collapse = ", "),
    call. = FALSE
  )
}
if (strict_packages && any(!tested_match)) {
  observed <- ifelse(
    tested_available,
    versions[names(tested_versions)],
    "not installed"
  )
  stop(
    "Strict validation requires ",
    paste(
      paste0(names(tested_versions), " ", tested_versions),
      collapse = " and "
    ),
    "; observed ",
    paste(paste0(names(observed), " ", observed), collapse = " and "),
    ".",
    call. = FALSE
  )
}

cat("R-code structure and 1,000-replicate configuration: PASS\n")
cat("V2.2 PDMI: B=", pdmi_config$B_sim,
    ", M=", pdmi_config$M_imp_pdmi, "\n", sep = "")
cat("Comparators: B=", comparator_config$B_sim,
    ", M_SMCFCS=", comparator_config$M_imp_smcfcs, "\n", sep = "")
cat("\nInstalled package check:\n")
print(package_table, row.names = FALSE)
if (any(!available) || any(!tested_match)) {
  cat(
    "\nRun scripts/install_dependencies.R before simulation. ",
    "Use --strict-packages to enforce the tested package versions.\n",
    sep = ""
  )
}
