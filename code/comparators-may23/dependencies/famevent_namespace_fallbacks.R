## ============================================================
## FamEvent namespace compatibility helpers.
##
## Some objects used by the related-work source files are internal
## FamEvent namespace objects and are not present in this local
## dependency dump. Import missing objects from an installed FamEvent
## package when possible, without modifying the original dependency files.
## ============================================================

.famevent_target_env <- environment()

.famevent_import_if_missing <- function(symbols, target_env = .famevent_target_env) {
  if (!requireNamespace("FamEvent", quietly = TRUE)) return(invisible(FALSE))
  ns <- asNamespace("FamEvent")
  for (sym in symbols) {
    if (!exists(sym, envir = target_env, inherits = TRUE) &&
        exists(sym, envir = ns, inherits = FALSE)) {
      assign(sym, get(sym, envir = ns, inherits = FALSE), envir = target_env)
    }
  }
  invisible(TRUE)
}

.famevent_import_if_missing(c(
  "cumhaz", "hazards", "gh", "penmodel", "loglik_frailty",
  "dlaplace", "laplace", "familyDesign", "fgeneZ", "fgeneZX",
  "Pgene", "surv.dist", "survp.dist", "inv.surv", "inv.survp",
  "inv2.surv", "parents.g", "kids.g", "simfam"
))

if (!exists("inv.survp", envir = .famevent_target_env, inherits = TRUE) &&
    exists("survp.dist", envir = .famevent_target_env, inherits = TRUE)) {
  inv.survp <- function(val, base.dist, parms, alpha) {
    out <- try(
      uniroot(
        survp.dist,
        lower = 0, upper = 100000,
        base.dist = base.dist,
        currentage = val[2],
        parms = parms,
        xbeta = val[1],
        alpha = alpha,
        res = val[3]
      )$root,
      silent = TRUE
    )
    if (is.null(attr(out, "class"))) return(out)
    print(c(parms, val))
  }
}

rm(.famevent_import_if_missing, .famevent_target_env)
