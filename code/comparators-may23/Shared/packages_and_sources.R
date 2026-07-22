## Author: Jiaqi Bi

## ============================================================
## Package and legacy-source loading for V2.
## Legacy files are loaded only for family skeleton/FamEvent support and
## kinship mechanics. V2 files override the missing-data and MLMI methods.
## ============================================================

v2_set_thread_env <- function() {
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    BLAS_NUM_THREADS = "1",
    LAPACK_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )
  invisible(TRUE)
}

v2_require_packages <- function(include_frailtypack = TRUE, include_smcfcs = FALSE) {
  required <- c("Matrix", "parallel", "survival", "kinship2", "MASS", "truncnorm")
  if (isTRUE(include_frailtypack)) required <- c(required, "frailtypack")
  if (isTRUE(include_smcfcs)) required <- c(required, "smcfcs")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Required package(s) not installed: ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages({
    library(Matrix)
    library(parallel)
    library(survival)
    library(kinship2)
    library(MASS)
    library(truncnorm)
    if (isTRUE(include_frailtypack)) library(frailtypack)
  })
  invisible(TRUE)
}

v2_source_if_exists <- function(path) {
  if (file.exists(path)) {
    source(path)
    return(TRUE)
  }
  FALSE
}

v2_load_family_sources <- function(code_root = v2_find_code_root(getwd())) {
  dep_dir <- file.path(code_root, "dependencies")
  rels <- c(
    "Delete Males.R",
    "FamEvent/R/cumhaz.R",
    "FamEvent/R/hazards.R",
    "FamEvent/R/gh.R",
    "FamEvent/R/penmodel.R",
    "FamEvent/R/loglik_frailty.R",
    "FamEvent/R/dlaplace.R",
    "FamEvent/R/laplace.R",
    "FamEvent/R/familyDesign.R",
    "FamEvent/R/fgeneZX.R",
    "FamEvent/R/Pgene.R",
    "FamEvent/R/surv.dist.R",
    "FamEvent/R/survp.dist.R",
    "FamEvent/R/inv.surv.R",
    "FamEvent/R/inv2.surv.R",
    "FamEvent/R/parents.g.R",
    "FamEvent/R/kids.g.R",
    "FamEvent/R/simfam.R",
    "famevent_namespace_fallbacks.R",
    "familyStructure_1to20.R"
  )
  loaded <- vapply(file.path(dep_dir, rels), v2_source_if_exists, logical(1))
  invisible(loaded)
}
