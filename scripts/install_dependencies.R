## Author: Jiaqi Bi

## Install the R packages needed by the continuous-covariate simulations.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(file_arg)) {
  sub("^--file=", "", file_arg[1])
} else {
  "scripts/install_dependencies.R"
}
script_file <- normalizePath(script_file, mustWork = TRUE)
repository_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)

local_library <- Sys.getenv(
  "REPRO_R_LIBRARY",
  file.path(repository_root, ".R-library")
)
dir.create(local_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(normalizePath(local_library, mustWork = TRUE), .libPaths())))

repos <- getOption("repos")
if (!length(repos) || identical(unname(repos[["CRAN"]]), "@CRAN@")) {
  repos <- c(CRAN = "https://cloud.r-project.org")
}

required_packages <- c(
  "Matrix", "survival", "kinship2", "MASS", "truncnorm",
  "frailtypack", "smcfcs"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages)) {
  message(
    "Installing packages into ", local_library, ": ",
    paste(missing_packages, collapse = ", ")
  )
  install.packages(
    missing_packages,
    lib = local_library,
    repos = repos,
    dependencies = c("Depends", "Imports", "LinkingTo")
  )
}

unavailable <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(unavailable)) {
  stop(
    "Package installation did not complete for: ",
    paste(unavailable, collapse = ", "),
    call. = FALSE
  )
}

tested_versions <- c(frailtypack = "2.13", smcfcs = "2.0.2")
observed_versions <- vapply(
  names(tested_versions),
  function(package) as.character(utils::packageVersion(package)),
  character(1)
)
version_mismatch <- observed_versions != tested_versions
if (any(version_mismatch)) {
  warning(
    "This code was validated with ",
    paste(
      paste0(names(tested_versions), " ", tested_versions),
      collapse = " and "
    ),
    ". Installed versions are ",
    paste(
      paste0(names(observed_versions), " ", observed_versions),
      collapse = " and "
    ),
    ".",
    call. = FALSE
  )
}

message("Dependency installation completed.")
message("R library: ", local_library)
print(sessionInfo())
