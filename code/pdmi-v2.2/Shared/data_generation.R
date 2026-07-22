## Author: Jiaqi Bi

## ============================================================
## Complete-data generator aligned to the V2.2 technical note.
## Complete families and covariates are generated first; missingness is
## applied only by Shared/missingness_and_cca.R.
## ============================================================

v22_make_kinship <- function(dat, id_col = "indID", father_col = "fatherID", mother_col = "motherID") {
  ids <- as.character(dat[[id_col]])
  Kraw <- 2 * kinship2::kinship(
    id = ids,
    dadid = as.character(dat[[father_col]]),
    momid = as.character(dat[[mother_col]])
  )
  K <- as.matrix(Kraw)
  K <- K[ids, ids, drop = FALSE]
  rownames(K) <- colnames(K) <- ids
  0.5 * (K + t(K))
}

v22_generate_mendelian_genotypes <- function(dat, pm = 0.02,
                                            id_col = "indID",
                                            father_col = "fatherID",
                                            mother_col = "motherID") {
  ids <- as.character(dat[[id_col]])
  father <- as.character(dat[[father_col]])
  mother <- as.character(dat[[mother_col]])
  father[father %in% c("", "0", "NA")] <- NA_character_
  mother[mother %in% c("", "0", "NA")] <- NA_character_
  fi <- match(father, ids)
  mi <- match(mother, ids)
  G <- rep(NA_integer_, length(ids))
  remaining <- seq_along(ids)
  guard <- 0L
  while (length(remaining)) {
    guard <- guard + 1L
    if (guard > length(ids) + 5L) {
      G[remaining] <- stats::rbinom(length(remaining), size = 2L, prob = pm)
      break
    }
    progressed <- FALSE
    for (loc in remaining) {
      has_parents <- !is.na(fi[loc]) && !is.na(mi[loc])
      parents_ready <- has_parents && !is.na(G[fi[loc]]) && !is.na(G[mi[loc]])
      if (!has_parents) {
        G[loc] <- stats::rbinom(1L, size = 2L, prob = pm)
        progressed <- TRUE
      } else if (parents_ready) {
        G[loc] <- stats::rbinom(1L, size = 1L, prob = G[fi[loc]] / 2) +
          stats::rbinom(1L, size = 1L, prob = G[mi[loc]] / 2)
        progressed <- TRUE
      }
    }
    remaining <- which(is.na(G))
    if (!progressed && length(remaining)) {
      G[remaining] <- stats::rbinom(length(remaining), size = 2L, prob = pm)
      break
    }
  }
  G
}

v22_generate_prs <- function(dat, K, sigma2_prs = 0.1) {
  out <- numeric(nrow(dat))
  blocks <- v22_family_blocks(dat)
  for (idx in blocks) {
    Ki <- as.matrix(K[idx, idx, drop = FALSE])
    out[idx] <- v22_rmvnorm_cov(rep(0, length(idx)), sigma2_prs * Ki)
  }
  out
}

v22_regenerate_event_times <- function(dat, K, omega, config = v22_default_config()) {
  th <- v22_theta_from_omega(omega)
  U <- numeric(nrow(dat))
  blocks <- v22_family_blocks(dat)
  for (idx in blocks) {
    Ki <- as.matrix(K[idx, idx, drop = FALSE])
    U[idx] <- v22_rmvnorm_cov(rep(0, length(idx)), th$sigma_u2 * Ki)
  }
  eta <- th$beta_c * as.numeric(dat$newx) + th$beta_b * as.numeric(dat$mgene) + U
  event_time <- th$lambda * (-log(stats::runif(nrow(dat))) / v22_safe_exp(eta))^(1 / th$rho)
  currentage <- pmin(pmax(as.numeric(dat$currentage), config$agemin), config$agemax)
  dat$u_true <- U
  dat$ageonset <- event_time
  dat$time <- pmin(currentage, event_time, config$agemax)
  dat$status <- as.integer(event_time <= currentage & event_time <= config$agemax)
  dat$t0 <- 0
  dat$fsize <- ave(dat$indID, dat$famID, FUN = length)
  dat
}

v22_prepare_pedigree_copy <- function(full_dat, analysis_dat) {
  ped <- full_dat
  for (nm in c("newx", "time", "status", "t0", "u_true", "ageonset", "fsize")) {
    if (!nm %in% names(ped)) ped[[nm]] <- NA_real_
  }
  ped$newx <- 0
  ped$time <- 0
  ped$status <- 0L
  ped$t0 <- 0
  ped$u_true <- 0
  ped$ageonset <- NA_real_
  ped$fsize <- ave(ped$indID, ped$famID, FUN = length)
  m <- match(as.character(analysis_dat$indID), as.character(ped$indID))
  for (nm in intersect(names(analysis_dat), names(ped))) {
    ped[[nm]][m] <- analysis_dat[[nm]]
  }
  ped
}

v22_generate_candidate_family <- function(fam_id, cumind, sigma_u2, config = v22_default_config()) {
  skel <- familyStructure(
    i = fam_id,
    cumind = cumind,
    m.carrier = 0,
    variation = config$variation,
    interaction = config$interaction,
    add.x = TRUE,
    x.dist = "normal",
    x.parms = c(0, sqrt(config$prs_sigma2), 1),
    depend = sigma_u2^(-1),
    base.dist = config$base.dist,
    frailty.dist = config$frailty.dist,
    base.parms = c(exp(-unname(config$omega_base["log.lambda"])),
                   exp(unname(config$omega_base["log.rho"]))),
    vbeta = c(0, unname(config$omega_base["beta_b"]),
              unname(config$omega_base["beta_c"])),
    allelefreq = config$pm,
    dominant.m = TRUE,
    dominant.s = TRUE,
    mrate = 0,
    probandage = config$probandage,
    agemin = config$agemin,
    agemax = config$agemax,
    align_female_only = config$align_female_only
  )
  full <- as.data.frame(skel)
  full$geno_A_count <- v22_generate_mendelian_genotypes(full, pm = config$pm)
  full$mgene <- as.integer(full$geno_A_count >= 1L)
  full$majorgene <- 3L - full$geno_A_count

  analysis <- if (isTRUE(config$analysis_female_only) && "gender" %in% names(full)) {
    full[as.numeric(full$gender) == 0, , drop = FALSE]
  } else {
    full
  }
  analysis <- analysis[order(as.numeric(analysis$indID)), , drop = FALSE]
  rownames(analysis) <- NULL
  K <- v22_make_kinship(analysis)
  analysis$newx <- v22_generate_prs(analysis, K, sigma2_prs = config$prs_sigma2)
  omega <- v22_actual_omega(sigma_u2, config)
  analysis <- v22_regenerate_event_times(analysis, K, omega, config)
  ped <- v22_prepare_pedigree_copy(full, analysis)
  list(dat = analysis, K = K, pedigree_dat = ped, full_dat = full)
}

v22_check_popplus_support <- function(dat) {
  blocks <- v22_family_blocks(dat)
  checks <- lapply(blocks, function(idx) {
    p <- which(as.numeric(dat$proband[idx]) == 1)
    if (length(p) != 1L) return(FALSE)
    j <- idx[p]
    isTRUE(as.numeric(dat$mgene[j]) == 1) && isTRUE(as.numeric(dat$status[j]) == 1)
  })
  all(unlist(checks))
}

v22_generate_complete_data <- function(replicate_id, sigma_u2, seed,
                                      config = v22_default_config()) {
  set.seed(seed)
  out <- vector("list", config$n_families)
  K_out <- vector("list", config$n_families)
  ped_out <- vector("list", config$n_families)
  full_out <- vector("list", config$n_families)
  fam_id <- 1L
  cumind <- 0L
  attempts <- 0L
  attempts_per_family <- config$selection_max_attempts_per_family %||% 5000L
  max_attempts <- max(1000L, as.integer(attempts_per_family) * config$n_families)
  while (fam_id <= config$n_families) {
    attempts <- attempts + 1L
    if (attempts > max_attempts) {
      stop("Exceeded max_attempts while generating selected pop+ families.")
    }
    cand <- tryCatch(v22_generate_candidate_family(fam_id, cumind, sigma_u2, config),
                     error = function(e) NULL)
    if (is.null(cand)) next
    prob <- which(as.numeric(cand$dat$proband) == 1)
    if (length(prob) != 1L) next
    keep <- isTRUE(as.numeric(cand$dat$mgene[prob]) == 1) &&
      isTRUE(as.numeric(cand$dat$status[prob]) == 1)
    if (!keep) next
    out[[fam_id]] <- cand$dat
    K_out[[fam_id]] <- cand$K
    ped_out[[fam_id]] <- cand$pedigree_dat
    full_out[[fam_id]] <- cand$full_dat
    cumind <- cumind + nrow(cand$full_dat)
    fam_id <- fam_id + 1L
  }
  dat <- do.call(rbind, out)
  rownames(dat) <- NULL
  K <- as.matrix(Matrix::bdiag(lapply(K_out, as.matrix)))
  rownames(K) <- colnames(K) <- as.character(dat$indID)
  pedigree_dat <- do.call(rbind, ped_out)
  rownames(pedigree_dat) <- NULL
  full_dat <- do.call(rbind, full_out)
  rownames(full_dat) <- NULL
  if (!v22_check_popplus_support(dat)) stop("Generated data violate pop+ support.")
  list(
    replicate_id = replicate_id,
    sigma_u2 = sigma_u2,
    dat = dat,
    K = v22_align_K(K, dat),
    pedigree_dat = pedigree_dat,
    full_dat = full_dat,
    generation_attempts = attempts
  )
}
