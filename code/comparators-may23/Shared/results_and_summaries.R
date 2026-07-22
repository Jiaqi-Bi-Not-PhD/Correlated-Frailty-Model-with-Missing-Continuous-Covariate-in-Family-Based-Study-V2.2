## Author: Jiaqi Bi

## ============================================================
## Result packaging and performance summaries.
## ============================================================

v2_pack_result <- function(metadata, fit = NULL, pooled = NULL, penetrance = NULL,
                           diagnostics = list(), convergence = TRUE,
                           failure_reason = NA_character_) {
  omega <- if (!is.null(pooled)) pooled$omega else if (!is.null(fit)) fit$omega else setNames(rep(NA_real_, 5), v2_omega_names())
  V <- if (!is.null(pooled)) pooled$vcov_omega else if (!is.null(fit)) fit$vcov_omega else matrix(NA_real_, 5, 5)
  dimnames(V) <- list(v2_omega_names(), v2_omega_names())
  list(
    metadata = metadata,
    convergence = isTRUE(convergence),
    failure_reason = failure_reason,
    omega = omega[v2_omega_names()],
    vcov_omega = V[v2_omega_names(), v2_omega_names(), drop = FALSE],
    penetrance = penetrance,
    diagnostics = diagnostics,
    pooled = pooled,
    fit = fit
  )
}

v2_result_path <- function(metadata, config = v2_default_config()) {
  rate <- if (is.na(metadata$target_missing_rate)) "nomiss" else sprintf("miss%02d", round(100 * metadata$target_missing_rate))
  file.path(
    config$results_root,
    metadata$run_label,
    metadata$missing_type,
    paste0("sigma_u2_", gsub("\\.", "p", as.character(metadata$sigma_u2))),
    rate,
    v2_method_tag(metadata$method, metadata$prior_version),
    sprintf("replicate_%04d.rds", metadata$replicate_id)
  )
}

v2_is_once_only_benchmark <- function(metadata) {
  method <- tolower(metadata$method %||% "")
  method %in% c("full-data", "cca", "mi-smcfcs")
}

v2_existing_result_path <- function(metadata, config = v2_default_config()) {
  path <- v2_result_path(metadata, config)
  if (isTRUE(config$skip_existing_results) &&
      file.exists(path) &&
      isTRUE(tryCatch({ readRDS(path); TRUE }, error = function(e) FALSE))) {
    return(path)
  }
  NA_character_
}

v2_existing_benchmark_path <- function(metadata, config = v2_default_config()) {
  path <- v2_result_path(metadata, config)
  if (isTRUE(config$skip_existing_benchmarks) &&
      v2_is_once_only_benchmark(metadata) &&
      file.exists(path)) {
    return(path)
  }
  NA_character_
}

v2_write_result <- function(result, config = v2_default_config()) {
  path <- v2_result_path(result$metadata, config)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, path)
  path
}

v2_summary_group_key <- function(...) {
  vals <- lapply(list(...), function(x) {
    x <- as.character(x)
    x[is.na(x)] <- "NA"
    x
  })
  do.call(paste, c(vals, sep = "\r"))
}

v2_scenario_label <- function(missing_type, sigma_u2, target_missing_rate, prior_version = NA_character_) {
  miss <- if (is.na(target_missing_rate)) "no missing" else paste0("target missing ", round(100 * target_missing_rate), "%")
  prior <- if (!is.na(prior_version) && nzchar(prior_version) && prior_version != "none") {
    paste0("; prior=", prior_version)
  } else ""
  paste0(missing_type, "; sigma_u2=", sigma_u2, "; ", miss, prior)
}

v2_rows_from_result <- function(result, true_omega = NULL) {
  md <- result$metadata
  true_omega <- true_omega %||% v2_actual_omega(md$sigma_u2)
  se <- sqrt(pmax(diag(result$vcov_omega), 0))
  data.frame(
    replicate_id = md$replicate_id,
    missing_type = md$missing_type,
    sigma_u2 = md$sigma_u2,
    target_missing_rate = md$target_missing_rate,
    method = md$method,
    prior_version = md$prior_version,
    M_imp = md$M_imp,
    parameter = v2_omega_names(),
    true = as.numeric(true_omega[v2_omega_names()]),
    estimate = as.numeric(result$omega[v2_omega_names()]),
    se = as.numeric(se[v2_omega_names()]),
    convergence = result$convergence,
    failure_reason = result$failure_reason,
    complete_data_seed = md$complete_data_seed,
    missing_mask_seed = md$missing_mask_seed,
    method_seed = md$method_seed,
    run_label = md$run_label,
    stringsAsFactors = FALSE
  )
}

v2_summarize_parameter_rows <- function(rows) {
  ok <- rows[is.finite(rows$estimate) & is.finite(rows$se) & rows$convergence, , drop = FALSE]
  if (!nrow(ok)) return(data.frame())
  split_key <- v2_summary_group_key(ok$missing_type, ok$sigma_u2, ok$target_missing_rate,
                                    ok$method, ok$prior_version, ok$parameter)
  pieces <- split(ok, split_key)
  do.call(rbind, lapply(pieces, function(d) {
    covered <- d$estimate - 1.96 * d$se <= d$true & d$estimate + 1.96 * d$se >= d$true
    emp_se <- stats::sd(d$estimate)
    model_se <- sqrt(mean(d$se^2))
    data.frame(
      scenario = v2_scenario_label(d$missing_type[1], d$sigma_u2[1],
                                   d$target_missing_rate[1], d$prior_version[1]),
      missing_type = d$missing_type[1],
      sigma_u2 = d$sigma_u2[1],
      target_missing_rate = d$target_missing_rate[1],
      method = d$method[1],
      prior_version = d$prior_version[1],
      parameter = d$parameter[1],
      true = d$true[1],
      n_success = nrow(d),
      bias = mean(d$estimate - d$true),
      rmse = sqrt(mean((d$estimate - d$true)^2)),
      model_se = model_se,
      empirical_se = emp_se,
      coverage = mean(covered),
      mcse_bias = emp_se / sqrt(nrow(d)),
      mcse_model_se = stats::sd(d$se) / sqrt(nrow(d)),
      mcse_empirical_se = if (nrow(d) > 1L) emp_se / sqrt(2 * (nrow(d) - 1L)) else NA_real_,
      mcse_coverage = sqrt(mean(covered) * (1 - mean(covered)) / nrow(d)),
      stringsAsFactors = FALSE
    )
  }))
}

v2_penetrance_rows_from_result <- function(result, config = v2_default_config()) {
  if (is.null(result$penetrance) || !nrow(result$penetrance)) return(data.frame())
  md <- result$metadata
  true_grid <- v2_penetrance_grid(v2_actual_omega(md$sigma_u2, config), config)
  key <- paste(true_grid$age, true_grid$prs, true_grid$gene)
  true_val <- true_grid$estimate[match(paste(result$penetrance$age, result$penetrance$prs, result$penetrance$gene), key)]
  data.frame(
    replicate_id = md$replicate_id,
    missing_type = md$missing_type,
    sigma_u2 = md$sigma_u2,
    target_missing_rate = md$target_missing_rate,
    method = md$method,
    prior_version = md$prior_version,
    M_imp = md$M_imp,
    age = result$penetrance$age,
    prs = result$penetrance$prs,
    gene = result$penetrance$gene,
    true = true_val,
    estimate = result$penetrance$estimate,
    se = result$penetrance$se,
    convergence = result$convergence,
    failure_reason = result$failure_reason,
    run_label = md$run_label,
    stringsAsFactors = FALSE
  )
}

v2_summarize_penetrance_rows <- function(rows) {
  ok <- rows[is.finite(rows$estimate) & is.finite(rows$se) & rows$convergence, , drop = FALSE]
  if (!nrow(ok)) return(data.frame())
  split_key <- v2_summary_group_key(ok$missing_type, ok$sigma_u2, ok$target_missing_rate,
                                    ok$method, ok$prior_version, ok$age, ok$prs, ok$gene)
  pieces <- split(ok, split_key)
  do.call(rbind, lapply(pieces, function(d) {
    covered <- d$estimate - 1.96 * d$se <= d$true & d$estimate + 1.96 * d$se >= d$true
    emp_se <- stats::sd(d$estimate)
    model_se <- sqrt(mean(d$se^2))
    data.frame(
      scenario = v2_scenario_label(d$missing_type[1], d$sigma_u2[1],
                                   d$target_missing_rate[1], d$prior_version[1]),
      missing_type = d$missing_type[1],
      sigma_u2 = d$sigma_u2[1],
      target_missing_rate = d$target_missing_rate[1],
      method = d$method[1],
      prior_version = d$prior_version[1],
      age = d$age[1],
      prs = d$prs[1],
      gene = d$gene[1],
      true = d$true[1],
      n_success = nrow(d),
      bias = mean(d$estimate - d$true),
      rmse = sqrt(mean((d$estimate - d$true)^2)),
      model_se = model_se,
      empirical_se = emp_se,
      coverage = mean(covered),
      mcse_bias = emp_se / sqrt(nrow(d)),
      mcse_model_se = stats::sd(d$se) / sqrt(nrow(d)),
      mcse_empirical_se = if (nrow(d) > 1L) emp_se / sqrt(2 * (nrow(d) - 1L)) else NA_real_,
      mcse_coverage = sqrt(mean(covered) * (1 - mean(covered)) / nrow(d)),
      stringsAsFactors = FALSE
    )
  }))
}

v2_manifest_parameter_table <- function(param_summary) {
  if (!nrow(param_summary)) return(data.frame())
  data.frame(
    Scenario = param_summary$scenario,
    Method = ifelse(param_summary$prior_version %in% c("none", NA_character_),
                    param_summary$method,
                    paste(param_summary$method, param_summary$prior_version)),
    Parameter = param_summary$parameter,
    True_value = param_summary$true,
    Bias = param_summary$bias,
    RMSE = param_summary$rmse,
    Model_SE = param_summary$model_se,
    Empirical_SE = param_summary$empirical_se,
    Coverage = param_summary$coverage,
    MCSE_Bias = param_summary$mcse_bias,
    MCSE_Model_SE = param_summary$mcse_model_se,
    MCSE_Empirical_SE = param_summary$mcse_empirical_se,
    MCSE_Coverage = param_summary$mcse_coverage,
    stringsAsFactors = FALSE
  )
}

v2_manifest_penetrance_table <- function(pen_summary) {
  if (!nrow(pen_summary)) return(data.frame())
  data.frame(
    Scenario = pen_summary$scenario,
    Age = pen_summary$age,
    PRS = pen_summary$prs,
    Major_gene_profile = pen_summary$gene,
    Method = ifelse(pen_summary$prior_version %in% c("none", NA_character_),
                    pen_summary$method,
                    paste(pen_summary$method, pen_summary$prior_version)),
    True_value = pen_summary$true,
    Bias = pen_summary$bias,
    RMSE = pen_summary$rmse,
    Model_SE = pen_summary$model_se,
    Empirical_SE = pen_summary$empirical_se,
    Coverage = pen_summary$coverage,
    MCSE_Bias = pen_summary$mcse_bias,
    MCSE_Model_SE = pen_summary$mcse_model_se,
    MCSE_Empirical_SE = pen_summary$mcse_empirical_se,
    MCSE_Coverage = pen_summary$mcse_coverage,
    stringsAsFactors = FALSE
  )
}

v2_flatten_diagnostics <- function(x, prefix = "diagnostic") {
  out <- list()
  walk <- function(y, nm) {
    if (is.list(y) && !is.data.frame(y)) {
      nms <- names(y) %||% paste0("item", seq_along(y))
      for (k in seq_along(y)) walk(y[[k]], paste(nm, nms[k], sep = "."))
    } else if (length(y) == 1L && (is.atomic(y) || is.null(y))) {
      out[[nm]] <<- if (is.null(y)) NA else y
    }
  }
  walk(x, prefix)
  out
}

v2_diagnostic_rows_from_results <- function(results) {
  rows <- lapply(results, function(result) {
    md <- result$metadata
    base <- list(
      replicate_id = md$replicate_id,
      missing_type = md$missing_type,
      sigma_u2 = md$sigma_u2,
      target_missing_rate = md$target_missing_rate,
      method = md$method,
      prior_version = md$prior_version,
      M_imp = md$M_imp,
      convergence = result$convergence,
      failure_reason = result$failure_reason,
      complete_data_seed = md$complete_data_seed,
      missing_mask_seed = md$missing_mask_seed,
      method_seed = md$method_seed,
      run_label = md$run_label
    )
    as.data.frame(c(base, v2_flatten_diagnostics(result$diagnostics)), stringsAsFactors = FALSE)
  })
  Reduce(function(a, b) {
    all_names <- sort(unique(c(names(a), names(b))))
    miss_a <- setdiff(all_names, names(a))
    miss_b <- setdiff(all_names, names(b))
    for (nm in miss_a) a[[nm]] <- NA
    for (nm in miss_b) b[[nm]] <- NA
    rbind(a[all_names], b[all_names])
  }, rows)
}

v2_summarize_diagnostics <- function(diagnostic_rows) {
  if (!nrow(diagnostic_rows)) return(data.frame())
  split_key <- v2_summary_group_key(diagnostic_rows$missing_type, diagnostic_rows$sigma_u2,
                                    diagnostic_rows$target_missing_rate, diagnostic_rows$method,
                                    diagnostic_rows$prior_version)
  pieces <- split(diagnostic_rows, split_key)
  do.call(rbind, lapply(pieces, function(d) {
    failed <- d$failure_reason[!d$convergence]
    failed <- failed[!is.na(failed) & nzchar(failed)]
    common <- if (length(failed)) {
      paste(names(sort(table(failed), decreasing = TRUE))[seq_len(min(5L, length(table(failed))))],
            collapse = " | ")
    } else ""
    data.frame(
      scenario = v2_scenario_label(d$missing_type[1], d$sigma_u2[1],
                                   d$target_missing_rate[1], d$prior_version[1]),
      missing_type = d$missing_type[1],
      sigma_u2 = d$sigma_u2[1],
      target_missing_rate = d$target_missing_rate[1],
      method = d$method[1],
      prior_version = d$prior_version[1],
      n_replicates = nrow(d),
      n_success = sum(d$convergence),
      failure_rate = mean(!d$convergence),
      common_failure_reasons = common,
      stringsAsFactors = FALSE
    )
  }))
}

v2_penetrance_curve_means <- function(pen_rows) {
  ok <- pen_rows[is.finite(pen_rows$estimate) & pen_rows$convergence, , drop = FALSE]
  if (!nrow(ok)) return(data.frame())
  split_key <- v2_summary_group_key(ok$missing_type, ok$sigma_u2, ok$target_missing_rate,
                                    ok$method, ok$prior_version, ok$age, ok$prs, ok$gene)
  pieces <- split(ok, split_key)
  do.call(rbind, lapply(pieces, function(d) {
    data.frame(
      missing_type = d$missing_type[1],
      sigma_u2 = d$sigma_u2[1],
      target_missing_rate = d$target_missing_rate[1],
      method = d$method[1],
      prior_version = d$prior_version[1],
      age = d$age[1],
      prs = d$prs[1],
      gene = d$gene[1],
      true = d$true[1],
      mean_estimate = mean(d$estimate),
      n_success = length(unique(d$replicate_id)),
      stringsAsFactors = FALSE
    )
  }))
}

v2_penetrance_integrated_accuracy <- function(pen_rows) {
  ok <- pen_rows[is.finite(pen_rows$estimate) & pen_rows$convergence, , drop = FALSE]
  if (!nrow(ok)) return(data.frame())
  split_key <- v2_summary_group_key(ok$missing_type, ok$sigma_u2, ok$target_missing_rate,
                                    ok$method, ok$prior_version, ok$prs, ok$gene)
  pieces <- split(ok, split_key)
  do.call(rbind, lapply(pieces, function(d) {
    age_means <- stats::aggregate(cbind(estimate, true) ~ age, data = d, FUN = mean)
    data.frame(
      scenario = v2_scenario_label(d$missing_type[1], d$sigma_u2[1],
                                   d$target_missing_rate[1], d$prior_version[1]),
      missing_type = d$missing_type[1],
      sigma_u2 = d$sigma_u2[1],
      target_missing_rate = d$target_missing_rate[1],
      method = d$method[1],
      prior_version = d$prior_version[1],
      prs = d$prs[1],
      gene = d$gene[1],
      IAB = mean(abs(age_means$estimate - age_means$true)),
      IRMSE = sqrt(mean((d$estimate - d$true)^2)),
      stringsAsFactors = FALSE
    )
  }))
}

v2_clean_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  gsub("_+", "_", x)
}

v2_penetrance_plot_groups <- function(curve_means) {
  if (!nrow(curve_means)) return(list())
  scenarios <- unique(curve_means[c("missing_type", "sigma_u2", "target_missing_rate")])
  groups <- list()
  for (i in seq_len(nrow(scenarios))) {
    sc <- scenarios[i, , drop = FALSE]
    d <- curve_means[
      curve_means$missing_type == sc$missing_type &
        curve_means$sigma_u2 == sc$sigma_u2 &
        ((is.na(curve_means$target_missing_rate) & is.na(sc$target_missing_rate)) |
           curve_means$target_missing_rate == sc$target_missing_rate),
      , drop = FALSE
    ]
    if (!nrow(d)) next
    priors <- unique(d$prior_version[d$method == "MI-Cong" & !is.na(d$prior_version)])
    priors <- priors[nzchar(priors) & priors != "none"]
    if (!length(priors)) priors <- unique(d$prior_version)
    for (prior in priors) {
      keep <- d$prior_version == prior |
        d$method %in% c("Full-data", "full-data", "CCA", "MI-SMCFCS") |
        d$method == "full-data"
      groups[[length(groups) + 1L]] <- list(
        missing_type = sc$missing_type,
        sigma_u2 = sc$sigma_u2,
        target_missing_rate = sc$target_missing_rate,
        prior_version = prior,
        data = d[keep, , drop = FALSE]
      )
    }
  }
  groups
}

v2_write_penetrance_curve_figures <- function(curve_means, out_dir, config = v2_default_config()) {
  if (!nrow(curve_means)) return(character())
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  groups <- v2_penetrance_plot_groups(curve_means)
  written <- character()
  colors <- c("full-data" = "#1F78B4", "Full-data" = "#1F78B4",
              "CCA" = "#D95F02", "MI-SMCFCS" = "#7570B3",
              "MI-Cong" = "#1B9E77")
  pchs <- c("full-data" = 16, "Full-data" = 16, "CCA" = 17, "MI-SMCFCS" = 15, "MI-Cong" = 18)
  for (grp in groups) {
    d <- grp$data
    if (!nrow(d)) next
    title <- v2_scenario_label(grp$missing_type, grp$sigma_u2,
                               grp$target_missing_rate, grp$prior_version)
    stem <- v2_clean_filename(paste("penetrance_curves", title, sep = "_"))
    for (device in c("pdf", "png")) {
      path <- file.path(out_dir, paste0(stem, ".", device))
      if (device == "pdf") {
        grDevices::pdf(path, width = 10, height = 7)
      } else {
        grDevices::png(path, width = 1600, height = 1100, res = 150)
      }
      op <- par(mfrow = c(length(config$penetrance_gene), length(config$penetrance_prs)),
                mar = c(4, 4, 2.6, 1), oma = c(0, 0, 3, 0))
      for (gene in config$penetrance_gene) {
        for (prs in config$penetrance_prs) {
          panel <- d[d$gene == gene & d$prs == prs, , drop = FALSE]
          true_panel <- panel[!duplicated(panel$age), c("age", "true")]
          true_panel <- true_panel[order(true_panel$age), , drop = FALSE]
          ylim <- range(c(panel$mean_estimate, true_panel$true), finite = TRUE)
          if (!all(is.finite(ylim))) ylim <- c(0, 1)
          ylim <- range(c(0, ylim, 1))
          plot(config$penetrance_ages, rep(NA_real_, length(config$penetrance_ages)),
               xlim = range(config$penetrance_ages), ylim = ylim,
               xlab = "Age", ylab = "Penetrance",
               main = paste0("PRS=", prs, ", gene=", gene))
          if (nrow(true_panel)) {
            lines(true_panel$age, true_panel$true, lwd = 2.5,
                  col = grDevices::adjustcolor("black", alpha.f = 0.5))
          }
          methods <- unique(panel$method)
          for (method in methods) {
            mdat <- panel[panel$method == method, , drop = FALSE]
            mdat <- mdat[order(mdat$age), , drop = FALSE]
            col <- unname(colors[method])
            if (is.na(col)) col <- "#444444"
            pch <- unname(pchs[method])
            if (is.na(pch)) pch <- 16
            points(mdat$age, mdat$mean_estimate, pch = pch,
                   col = grDevices::adjustcolor(col, alpha.f = 0.35))
            lines(mdat$age, mdat$mean_estimate, lty = 2, lwd = 1.5,
                  col = grDevices::adjustcolor(col, alpha.f = 0.35))
          }
          legend_cols <- unname(colors[methods])
          legend_cols[is.na(legend_cols)] <- "#444444"
          legend_pchs <- unname(pchs[methods])
          legend_pchs[is.na(legend_pchs)] <- 16
          legend("topleft", bty = "n", cex = 0.75,
                 legend = c("Truth", methods),
                 col = c(grDevices::adjustcolor("black", alpha.f = 0.5),
                         grDevices::adjustcolor(legend_cols, alpha.f = 0.35)),
                 lty = c(1, rep(2, length(methods))),
                 pch = c(NA, legend_pchs))
        }
      }
      mtext(title, outer = TRUE, cex = 1.1, font = 2)
      par(op)
      grDevices::dev.off()
      written <- c(written, path)
    }
  }
  written
}
