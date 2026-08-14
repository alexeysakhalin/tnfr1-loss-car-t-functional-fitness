# Figure 5D: prespecified bulk T-cell expression score and OS in CheckMate.
# Primary population: NIVOLUMAB only (n = 181). Inference is stratified by
# trial (CM-009/CM-010/CM-025). OS_CNSR is used unchanged: 1=event, 0=censored.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(survival)
  library(survminer)
  library(broom)
})

source(file.path("R", "00_bioinfo_helpers.R"))
source(file.path("R", "00_load_analysis_tables.R"))
source(file.path("R", "plot_style.R"))

out_dir <- "figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tcell_genes <- TCELL_SCORE_GENES
analysis_long <- load_analysis_tables()
assert_columns(analysis_long, c("gene_symbol", "expr_value"), "analysis_long")

dat <- build_nivolumab_expression(analysis_long, tcell_genes)
if (nrow(dat) != 181L || anyDuplicated(dat$sample_id)) {
  stop("Figure 5D requires 181 unique nivolumab RNA samples.")
}
present <- intersect(tcell_genes, names(dat))
if (!setequal(present, tcell_genes)) {
  stop("Figure 5D requires all five prespecified T-cell score genes.")
}
if (any(!is.finite(as.matrix(dat[present])))) {
  stop("Figure 5D requires complete expression for all five score genes.")
}
dat <- dat |>
  mutate(across(all_of(present), safe_z)) |>
  mutate(Tcell_score = rowMeans(across(all_of(present)), na.rm = FALSE))

fit <- fit_trial_stratified_survival(dat, "Tcell_score", "OS", "Tcell")
write.csv(fit$result,
          file.path(out_dir, "Figure_5D_nivolumab_trial_stratified_results.csv"),
          row.names = FALSE)

r <- fit$result[1, ]
subtitle <- sprintf(
  "Nivolumab only; trial-stratified Cox HR %.2f (95%% CI %.2f-%.2f), p=%.3f; descriptive log-rank p=%.3f",
  r$HR_high_vs_low, r$CI_low, r$CI_high, r$cox_p, r$logrank_p
)
p <- survminer::ggsurvplot(
  fit$km_fit, data = fit$data, risk.table = TRUE, censor = TRUE,
  palette = c("grey65", "#D62728"),
  legend.title = "T-cell score", legend.labs = c("Low", "High"),
  title = "Overall survival by bulk T-cell expression score",
  subtitle = subtitle, xlab = "Time (months)", ylab = "Overall survival probability"
)
save_ggsurvplot_png(
  p,
  file.path(out_dir, "Figure_5D_KM_ccRCC_Tcell_score.png")
)

qc <- tibble(
  item = c("analysis_arm", "n", "OS_events", "event_coding", "trial_adjustment",
           "Tcell_genes_present"),
  value = c("NIVOLUMAB", nrow(fit$data), sum(fit$data$event == 1),
            "1=event;0=censored", "strata(trial_id)", paste(present, collapse = ";"))
)
write.csv(qc, file.path(out_dir, "Figure_5D_QC.csv"), row.names = FALSE)
