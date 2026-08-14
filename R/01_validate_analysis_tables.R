# Validate the compact, analysis-ready tables before reproducing the figures.
# The source downloads and deterministic extraction steps are documented in
# the repository data manifest; no database build is required.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

source(file.path("R", "00_load_analysis_tables.R"))

analysis_long <- load_analysis_tables()

expected_cohorts <- c(
  "IMvigor210_BLCA", "CHECKMATE_CCRCC",
  "SU2C_MARK_NSCLC", "LIU2019_MELANOMA"
)
missing_cohorts <- setdiff(expected_cohorts, unique(analysis_long$cohort_id))
if (length(missing_cohorts)) {
  stop("Missing analysis cohorts: ", paste(missing_cohorts, collapse = ", "))
}

qc <- analysis_long |>
  group_by(.data$cohort_id, .data$expression_unit) |>
  summarise(
    n_samples = n_distinct(.data$sample_id),
    n_genes = n_distinct(.data$gene_symbol),
    n_expression_rows = n(),
    rank_percentile_missing = sum(is.na(.data$rank_percentile)),
    .groups = "drop"
  ) |>
  arrange(match(.data$cohort_id, expected_cohorts))

dir.create("results", recursive = TRUE, showWarnings = FALSE)
write_tsv(qc, file.path("results", "analysis_table_qc.tsv"))
print(qc)
