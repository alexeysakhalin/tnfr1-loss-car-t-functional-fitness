# Figure 5C: cross-cohort bulk T-cell expression score.
#
# rank_percentile is calculated during cohort-specific data preparation after
# mapping and duplicate-feature aggregation, against every mapped gene in the
# corresponding sample. It must not be recalculated from the selected-gene
# table because that table contains only analysis genes.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

source(file.path("R", "00_load_analysis_tables.R"))
source(file.path("R", "00_bioinfo_helpers.R"))
source(file.path("R", "plot_style.R"))

figure_dir <- "figures"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
analysis_long <- load_analysis_tables()

tcell_genes <- TCELL_SCORE_GENES
cohort_order <- c(
  "IMvigor210_BLCA", "CHECKMATE_CCRCC",
  "SU2C_MARK_NSCLC", "LIU2019_MELANOMA"
)
cohort_labels <- c(
  IMvigor210_BLCA = "Bladder",
  CHECKMATE_CCRCC = "ccRCC",
  SU2C_MARK_NSCLC = "Lung",
  LIU2019_MELANOMA = "Melanoma"
)
cohort_colors <- c(
  Bladder = "#2C7BB6", ccRCC = "#D62728",
  Lung = "#F28E2B", Melanoma = "#4CAF50"
)

score_long <- analysis_long |>
  filter(
    .data$cohort_id %in% cohort_order,
    .data$gene_symbol %in% tcell_genes,
    .data$cohort_id != "CHECKMATE_CCRCC" |
      toupper(trimws(as.character(.data$treatment_arm))) == "NIVOLUMAB"
  )
missing_cohorts <- setdiff(cohort_order, unique(score_long$cohort_id))
if (length(missing_cohorts)) {
  stop(
    "Figure 5C is missing cohorts: ", paste(missing_cohorts, collapse = ", ")
  )
}

eligible_samples <- analysis_long |>
  filter(
    .data$cohort_id %in% cohort_order,
    .data$cohort_id != "CHECKMATE_CCRCC" |
      toupper(trimws(as.character(.data$treatment_arm))) == "NIVOLUMAB"
  ) |>
  distinct(.data$cohort_id, .data$sample_id)
coverage <- score_long |>
  group_by(.data$cohort_id, .data$sample_id) |>
  summarise(
    n_genes = n_distinct(.data$gene_symbol),
    n_finite_ranks = sum(is.finite(.data$rank_percentile)),
    .groups = "drop"
  )
missing_samples <- eligible_samples |>
  anti_join(coverage, by = c("cohort_id", "sample_id"))
if (nrow(missing_samples) ||
    any(coverage$n_genes != length(tcell_genes)) ||
    any(coverage$n_finite_ranks != length(tcell_genes))) {
  stop(
    "Figure 5C requires all five precomputed marker-gene ranks for every sample."
  )
}

plot_data <- score_long |>
  group_by(.data$cohort_id, .data$sample_id) |>
  summarise(
    Tcell_rank_score = mean(.data$rank_percentile),
    n_genes_used = n_distinct(.data$gene_symbol),
    .groups = "drop"
  ) |>
  mutate(
    cohort_label = factor(
      unname(cohort_labels[.data$cohort_id]),
      levels = unname(cohort_labels[cohort_order])
    )
  )

checkmate_n <- plot_data |>
  filter(.data$cohort_id == "CHECKMATE_CCRCC") |>
  summarise(n = n_distinct(.data$sample_id)) |>
  pull(.data$n)
if (!identical(checkmate_n, 181L)) {
  stop("Figure 5C requires 181 nivolumab-treated CheckMate samples.")
}

statistics <- plot_data |>
  group_by(.data$cohort_label) |>
  summarise(
    n = n(),
    median_score = median(.data$Tcell_rank_score),
    iqr = IQR(.data$Tcell_rank_score),
    maximum = max(.data$Tcell_rank_score),
    .groups = "drop"
  )
score_range <- diff(range(plot_data$Tcell_rank_score))
if (!is.finite(score_range) || score_range == 0) {
  stop("Figure 5C scores have no finite between-sample range.")
}
statistics <- statistics |>
  mutate(
    label = paste0(
      "median = ", sprintf("%.2f", .data$median_score), "\n",
      "IQR = ", sprintf("%.2f", .data$iqr), "\n",
      "n = ", .data$n
    ),
    label_y = .data$maximum + 0.16 * score_range
  )
note_y <- min(plot_data$Tcell_rank_score) - 0.06 * score_range

plot_object <- ggplot(
  plot_data,
  aes(.data$cohort_label, .data$Tcell_rank_score, fill = .data$cohort_label)
) +
  geom_violin(
    alpha = 0.6, width = 0.55, trim = FALSE,
    scale = "width", color = NA
  ) +
  geom_boxplot(width = 0.12, fill = "white", color = "black", linewidth = 1.1) +
  geom_text(
    data = statistics,
    aes(x = .data$cohort_label, y = .data$label_y, label = .data$label),
    inherit.aes = FALSE, size = 4.5, fontface = "bold"
  ) +
  annotate(
    "text", x = 0.7, y = note_y,
    label = paste0(
      "Mean of five within-sample transcriptome-wide percentile ranks;\n",
      "IQR describes within-cohort dispersion"
    ),
    hjust = 0, size = 3.8
  ) +
  scale_fill_manual(values = cohort_colors) +
  scale_y_continuous(
    limits = c(
      min(plot_data$Tcell_rank_score) - 0.15 * score_range,
      max(statistics$label_y) + 0.12 * score_range
    ),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    title = "Cross-cohort comparison of a bulk T-cell expression score",
    x = NULL, y = "Rank-based T-cell expression score"
  ) +
  coord_cartesian(clip = "off") +
  theme_pub() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 15, hjust = 1),
    plot.margin = margin(30, 20, 20, 20)
  )

ggsave(
  file.path(figure_dir, "Figure_5C_bulk_Tcell_expression_score_rank_based.png"),
  plot_object, width = 11, height = 6.5, dpi = 600, bg = "white"
)
