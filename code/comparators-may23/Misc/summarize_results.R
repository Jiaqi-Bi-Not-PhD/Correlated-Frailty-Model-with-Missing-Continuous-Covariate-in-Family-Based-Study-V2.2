## Author: Jiaqi Bi

## Aggregate raw V2 replicate RDS files into parameter, penetrance,
## diagnostic CSVs, and penetrance curve figures required by the
## simulation section of the technical note.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "Misc/summarize_results.R"
script_file <- gsub("~+~", " ", script_file, fixed = TRUE)
code_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = FALSE)
source(file.path(code_root, "Shared", "source_all.R"))

config <- v2_default_config()
config$results_root <- Sys.getenv("SIM_RESULTS_ROOT", config$results_root)
run_label <- Sys.getenv("SIM_RUN_LABEL", "")
root <- if (nzchar(run_label)) file.path(config$results_root, run_label) else config$results_root
files <- list.files(root, pattern = "replicate_[0-9]+\\.rds$", recursive = TRUE, full.names = TRUE)
if (!length(files)) stop("No replicate RDS files found under ", root)

res <- lapply(files, readRDS)
v2_bind_rows <- function(x) {
  x <- x[vapply(x, nrow, integer(1)) > 0L]
  if (!length(x)) return(data.frame())
  do.call(rbind, x)
}

param_rows <- v2_bind_rows(lapply(res, v2_rows_from_result))
pen_rows <- v2_bind_rows(lapply(res, v2_penetrance_rows_from_result, config = config))
diagnostic_rows <- v2_diagnostic_rows_from_results(res)
param_summary <- v2_summarize_parameter_rows(param_rows)
pen_summary <- v2_summarize_penetrance_rows(pen_rows)
diagnostic_summary <- v2_summarize_diagnostics(diagnostic_rows)
curve_means <- v2_penetrance_curve_means(pen_rows)
integrated_accuracy <- v2_penetrance_integrated_accuracy(pen_rows)

out_dir <- file.path(root, "summary_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(param_rows, file.path(out_dir, "parameter_estimates_per_replicate.csv"), row.names = FALSE)
utils::write.csv(param_summary, file.path(out_dir, "parameter_performance_summary.csv"), row.names = FALSE)
utils::write.csv(v2_manifest_parameter_table(param_summary), file.path(out_dir, "parameter_manuscript_table.csv"), row.names = FALSE)
utils::write.csv(pen_rows, file.path(out_dir, "penetrance_estimates_per_replicate.csv"), row.names = FALSE)
utils::write.csv(pen_summary, file.path(out_dir, "penetrance_performance_summary.csv"), row.names = FALSE)
utils::write.csv(v2_manifest_penetrance_table(pen_summary), file.path(out_dir, "penetrance_manuscript_table.csv"), row.names = FALSE)
utils::write.csv(diagnostic_rows, file.path(out_dir, "diagnostics_per_replicate.csv"), row.names = FALSE)
utils::write.csv(diagnostic_summary, file.path(out_dir, "convergence_diagnostic_summary.csv"), row.names = FALSE)
utils::write.csv(curve_means, file.path(out_dir, "penetrance_curve_means.csv"), row.names = FALSE)
utils::write.csv(integrated_accuracy, file.path(out_dir, "penetrance_integrated_accuracy.csv"), row.names = FALSE)

fig_dir <- file.path(root, "figures", "penetrance_curves")
figures <- v2_write_penetrance_curve_figures(curve_means, fig_dir, config)
utils::write.csv(data.frame(path = figures), file.path(out_dir, "penetrance_curve_figure_manifest.csv"), row.names = FALSE)
message("Wrote summaries to ", out_dir)
message("Wrote penetrance figures to ", fig_dir)
