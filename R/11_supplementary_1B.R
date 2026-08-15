# Supplementary Figure S1B from the tracked DepMap derivative.
#
# scripts/prepare_depmap_s1b.py validates and reduces the large third-party
# sources. This renderer uses only packages already locked for figure builds
# and cannot inject a DepMap release label while source pairing is unverified.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

input_file <- file.path(
  "data", "analysis", "depmap_s1b_eligible_models.tsv.gz"
)
qc_file <- file.path(
  "data", "analysis", "depmap_s1b_preparation_qc.json"
)
source_provenance_file <- file.path(
  "data", "analysis", "depmap_s1b_source_provenance.json"
)
summary_file <- file.path(
  "reference_results", "depmap_s1b_statistics.csv"
)
out_dir <- Sys.getenv(
  "DEPMAP_OUTPUT_DIR",
  unset = file.path("results", "supplementary_S1B")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expected_sizes <- c(
  input = 47514,
  qc = 1583,
  source_provenance = 1895,
  summary = 305
)
tracked_files <- c(
  input = input_file,
  qc = qc_file,
  source_provenance = source_provenance_file,
  summary = summary_file
)
missing_files <- tracked_files[!file.exists(tracked_files)]
if (length(missing_files)) {
  stop("Missing tracked S1B input: ", unname(missing_files[[1]]))
}
observed_sizes <- vapply(
  tracked_files,
  function(path) unname(file.info(path)$size),
  numeric(1)
)
if (!identical(as.numeric(observed_sizes), as.numeric(expected_sizes))) {
  stop("A tracked S1B input has an unexpected byte size.")
}

# This state is intentionally not configurable at render time. Source hashes,
# release status and same-release=null are checked by the Python preparation
# contract and repository tests. Confirmation requires a reviewed data update.
release_pair_status <- "unverified"
expression_release <- "DepMap Public 25Q2"
model_release_identity_status <- "unverified"

expected_columns <- c(
  "ProfileID",
  "ModelID",
  "OncotreePrimaryDisease",
  "RIPK3_log2_TPM_plus_1",
  "NLRP3_log2_TPM_plus_1",
  "RIPK3_below_threshold",
  "NLRP3_below_threshold",
  "threshold_category"
)
analysis <- fread(input_file, showProgress = FALSE)
if (!identical(names(analysis), expected_columns)) {
  stop("The tracked S1B derivative has an unexpected column schema.")
}

parse_flag <- function(x, field) {
  if (is.logical(x)) {
    if (anyNA(x)) stop(field, " contains a missing value.")
    return(x)
  }
  normalized <- toupper(trimws(as.character(x)))
  if (!all(normalized %in% c("TRUE", "FALSE"))) {
    stop(field, " must contain only TRUE or FALSE.")
  }
  normalized == "TRUE"
}

analysis[, ProfileID := trimws(as.character(ProfileID))]
analysis[, ModelID := trimws(as.character(ModelID))]
analysis[, OncotreePrimaryDisease := trimws(
  as.character(OncotreePrimaryDisease)
)]
analysis[, RIPK3 := suppressWarnings(as.numeric(RIPK3_log2_TPM_plus_1))]
analysis[, NLRP3 := suppressWarnings(as.numeric(NLRP3_log2_TPM_plus_1))]
analysis[, RIPK3_below_threshold := parse_flag(
  RIPK3_below_threshold, "RIPK3_below_threshold"
)]
analysis[, NLRP3_below_threshold := parse_flag(
  NLRP3_below_threshold, "NLRP3_below_threshold"
)]
analysis[, threshold_category := trimws(as.character(threshold_category))]

if (nrow(analysis) != 1591L ||
    any(analysis$ProfileID == "") || anyDuplicated(analysis$ProfileID) ||
    any(analysis$ModelID == "") || anyDuplicated(analysis$ModelID) ||
    !all(grepl("^ACH-[0-9]{6}$", analysis$ModelID)) ||
    any(analysis$OncotreePrimaryDisease == "") ||
    any(toupper(analysis$OncotreePrimaryDisease) == "NON-CANCEROUS") ||
    !all(is.finite(analysis$RIPK3)) || !all(is.finite(analysis$NLRP3))) {
  stop("The tracked S1B derivative fails its row, identifier or value contract.")
}

cutoff <- 0.5
if (!identical(analysis$RIPK3_below_threshold, analysis$RIPK3 < cutoff) ||
    !identical(analysis$NLRP3_below_threshold, analysis$NLRP3 < cutoff)) {
  stop("Stored threshold flags do not agree with the expression values.")
}

expected_category <- ifelse(
  analysis$RIPK3_below_threshold & analysis$NLRP3_below_threshold,
  "Both below threshold",
  ifelse(
    analysis$RIPK3_below_threshold,
    "RIPK3 below threshold only",
    ifelse(
      analysis$NLRP3_below_threshold,
      "NLRP3 below threshold only",
      "Neither below threshold"
    )
  )
)
if (!identical(analysis$threshold_category, expected_category)) {
  stop("Stored threshold categories do not agree with the threshold flags.")
}

expected_metric_counts <- c(
  RIPK3_below_threshold = 1003L,
  NLRP3_below_threshold = 1172L,
  both_below_threshold = 749L,
  RIPK3_below_threshold_only = 254L,
  NLRP3_below_threshold_only = 423L,
  neither_below_threshold = 165L
)
observed_metric_counts <- c(
  RIPK3_below_threshold = sum(analysis$RIPK3_below_threshold),
  NLRP3_below_threshold = sum(analysis$NLRP3_below_threshold),
  both_below_threshold = sum(
    analysis$RIPK3_below_threshold & analysis$NLRP3_below_threshold
  ),
  RIPK3_below_threshold_only = sum(
    analysis$RIPK3_below_threshold & !analysis$NLRP3_below_threshold
  ),
  NLRP3_below_threshold_only = sum(
    !analysis$RIPK3_below_threshold & analysis$NLRP3_below_threshold
  ),
  neither_below_threshold = sum(
    !analysis$RIPK3_below_threshold & !analysis$NLRP3_below_threshold
  )
)
if (!identical(
  as.integer(observed_metric_counts), as.integer(expected_metric_counts)
)) {
  stop("Observed threshold counts differ from the locked S1B contract.")
}

summary_table <- fread(summary_file, showProgress = FALSE)
if (nrow(summary_table) != 6L ||
    any(summary_table$denominator != 1591L) ||
    any(summary_table$cutoff_log2_tpm_plus_1 != cutoff) ||
    !identical(
      as.integer(summary_table$n), as.integer(expected_metric_counts)
    ) ||
    any(abs(
      summary_table$percent - round(100 * summary_table$n / 1591, 1)
    ) > 1e-12)) {
  stop("The tracked statistics table differs from the locked S1B contract.")
}

category_levels <- c(
  "Both below threshold",
  "RIPK3 below threshold only",
  "NLRP3 below threshold only",
  "Neither below threshold"
)
category_labels <- c(
  "Both below threshold" = "Both below: n = 749 (47.1%)",
  "RIPK3 below threshold only" = "RIPK3 below only: n = 254 (16.0%)",
  "NLRP3 below threshold only" = "NLRP3 below only: n = 423 (26.6%)",
  "Neither below threshold" = "Neither below: n = 165 (10.4%)"
)
category_colors <- setNames(
  c("#4C78A8", "#F58518", "#54A24B", "#B279A2"),
  unname(category_labels[category_levels])
)
analysis[, plot_category := factor(
  unname(category_labels[threshold_category]),
  levels = unname(category_labels[category_levels])
)]

x_upper <- max(analysis$RIPK3) * 1.02
y_upper <- max(analysis$NLRP3) * 1.02
panel <- ggplot(
  analysis,
  aes(x = RIPK3, y = NLRP3, color = plot_category)
) +
  geom_point(size = 1.5, alpha = 0.58) +
  geom_vline(xintercept = cutoff, linetype = "dashed", color = "grey25") +
  geom_hline(yintercept = cutoff, linetype = "dashed", color = "grey25") +
  scale_color_manual(values = category_colors, drop = FALSE) +
  coord_cartesian(
    xlim = c(0, x_upper), ylim = c(0, y_upper),
    expand = FALSE, clip = "off"
  ) +
  guides(color = guide_legend(
    nrow = 2, byrow = TRUE,
    override.aes = list(size = 3, alpha = 1)
  )) +
  labs(
    title = "RIPK3 and NLRP3 expression in DepMap cell-line models",
    subtitle = paste0(
      "Eligible models (n = 1,591); dashed lines, 0.5 log2(TPM+1)"
    ),
    x = "RIPK3 expression, log2(TPM+1)",
    y = "NLRP3 expression, log2(TPM+1)",
    color = NULL,
    caption = paste0(
      "Eligibility: ModelType = Cell Line; OncoTree primary disease ",
      "not annotated as Non-Cancerous."
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 8.5),
    legend.key.width = grid::unit(0.9, "lines"),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10.5),
    plot.caption = element_text(hjust = 0, size = 8.5, color = "grey30"),
    plot.margin = margin(8, 12, 8, 8)
  )

ggsave(
  file.path(out_dir, "Supplementary_Figure_S1B.png"),
  panel, width = 7.1, height = 6.4, dpi = 600, bg = "white"
)
ggsave(
  file.path(out_dir, "Supplementary_Figure_S1B.tiff"),
  panel, width = 7.1, height = 6.4, dpi = 600,
  compression = "lzw", bg = "white"
)

runtime_provenance <- data.table(
  field = c(
    "release_pair_status",
    "expression_release",
    "model_release_identity_status",
    "derived_file",
    "population",
    "n_models",
    "cutoff_log2_tpm_plus_1",
    names(observed_metric_counts),
    "tissue_origin_filter_applied"
  ),
  value = c(
    release_pair_status,
    expression_release,
    model_release_identity_status,
    basename(input_file),
    paste0(
      "DepMap cell-line models with a non-missing OncoTree primary-disease ",
      "label other than Non-Cancerous"
    ),
    as.character(nrow(analysis)),
    as.character(cutoff),
    as.character(observed_metric_counts),
    "FALSE"
  )
)
fwrite(
  runtime_provenance,
  file.path(out_dir, "Supplementary_Figure_S1B_runtime_provenance.tsv"),
  sep = "\t", quote = TRUE
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
