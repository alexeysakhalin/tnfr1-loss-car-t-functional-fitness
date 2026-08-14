# Supplementary Figures S6A and S6C: cross-cohort receptor expression and an
# exploratory view of the visually defined high-TNFRSF1A peak in melanoma.
# All inputs are compact analysis tables with explicit cohort-specific units.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

source(file.path("R", "00_load_analysis_tables.R"))
source(file.path("R", "plot_style.R"))

figure_dir <- "figures"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
analysis_long <- load_analysis_tables()

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

primary_rows <- analysis_long |>
  filter(
    .data$cohort_id %in% cohort_order,
    .data$cohort_id != "CHECKMATE_CCRCC" |
      toupper(trimws(as.character(.data$treatment_arm))) == "NIVOLUMAB"
  )

plot_receptor_distribution <- function(gene, filename) {
  plot_data <- primary_rows |>
    filter(.data$gene_symbol == gene, is.finite(.data$expr_value)) |>
    mutate(
      cohort_label = factor(
        unname(cohort_labels[.data$cohort_id]),
        levels = unname(cohort_labels[cohort_order])
      )
    )

  missing <- setdiff(cohort_order, unique(plot_data$cohort_id))
  if (length(missing)) {
    stop(gene, " is missing from cohorts: ", paste(missing, collapse = ", "))
  }
  checkmate_n <- plot_data |>
    filter(.data$cohort_id == "CHECKMATE_CCRCC") |>
    summarise(n = n_distinct(.data$sample_id)) |>
    pull(.data$n)
  if (!identical(checkmate_n, 181L)) {
    stop(gene, " panel requires 181 nivolumab-treated CheckMate samples.")
  }

  units <- plot_data |>
    distinct(.data$cohort_label, .data$expression_unit) |>
    count(.data$cohort_label, name = "n_units")
  if (any(units$n_units != 1L)) {
    stop(gene, " has more than one expression unit within a cohort.")
  }

  medians <- plot_data |>
    group_by(.data$cohort_label) |>
    summarise(median_expression = median(.data$expr_value), .groups = "drop")
  unit_labels <- plot_data |>
    distinct(.data$cohort_label, .data$expression_unit) |>
    mutate(label = gsub("_", " ", .data$expression_unit, fixed = TRUE))

  plot_object <- ggplot(
    plot_data,
    aes(.data$expr_value, fill = .data$cohort_label, color = .data$cohort_label)
  ) +
    geom_histogram(
      aes(y = after_stat(density)), bins = 30, alpha = 0.25,
      position = "identity", color = NA
    ) +
    geom_density(linewidth = 1.4, alpha = 0.85) +
    geom_vline(
      data = medians,
      aes(xintercept = .data$median_expression, color = .data$cohort_label),
      linetype = "dashed", linewidth = 1, show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = unit_labels,
      aes(x = -Inf, y = -Inf, label = .data$label),
      inherit.aes = FALSE, hjust = -0.05, vjust = -0.2,
      size = 4, fontface = "bold"
    ) +
    facet_wrap(~ cohort_label, scales = "free", nrow = 1) +
    scale_fill_manual(values = cohort_colors) +
    scale_color_manual(values = cohort_colors) +
    labs(
      title = paste(gene, "expression across cohorts"),
      x = "Expression (cohort-specific units)", y = "Density"
    ) +
    coord_cartesian(clip = "off") +
    theme_pub() +
    theme(legend.position = "none")

  ggsave(
    file.path(figure_dir, filename), plot_object,
    width = 14, height = 4.7, dpi = 600, bg = "white"
  )
}

plot_receptor_distribution(
  "TNFRSF1A", "Supplementary_Figure_S6A_TNFRSF1A_expression.png"
)
plot_receptor_distribution(
  "IFNGR1", "Supplementary_Figure_S6A_IFNGR1_expression.png"
)

# The cutoff below describes a visually selected peak in the original Liu
# normalized-expression scale. It is used only for exploratory visualization; no
# melanoma samples are removed from another analysis.
tnfrsf1a_cutoff <- 200
liu_required_genes <- c("TNFRSF1A", "IFNG", "CD8A", "ICAM1", "IRF1")

liu_long <- primary_rows |>
  filter(
    .data$cohort_id == "LIU2019_MELANOMA",
    .data$gene_symbol %in% liu_required_genes
  )

liu_units <- unique(liu_long$expression_unit)
liu_unit_key <- if (length(liu_units) == 1L) {
  gsub("[^a-z0-9]", "", tolower(liu_units[[1]]))
} else {
  NA_character_
}
if (length(liu_units) != 1L ||
    !(liu_unit_key %in% c("tpm", "normalizedcounts", "normalizedexpression"))) {
  stop(
    "The melanoma peak cutoff requires the Liu normalized-count expression scale."
  )
}

melanoma <- liu_long |>
  select(
    .data$sample_id, .data$gene_symbol, .data$expr_value, .data$tumor_purity
  ) |>
  pivot_wider(names_from = .data$gene_symbol, values_from = .data$expr_value)

missing_genes <- setdiff(liu_required_genes, names(melanoma))
if (length(missing_genes)) {
  stop("Missing melanoma genes: ", paste(missing_genes, collapse = ", "))
}
if (anyDuplicated(melanoma$sample_id)) {
  stop("The melanoma analysis table contains duplicate samples.")
}

melanoma <- melanoma |>
  filter(is.finite(.data$TNFRSF1A)) |>
  mutate(
    right_peak_flag = .data$TNFRSF1A > tnfrsf1a_cutoff,
    peak_group = factor(
      if_else(.data$right_peak_flag, "Right peak", "Left peak"),
      levels = c("Left peak", "Right peak")
    )
  )
if (n_distinct(melanoma$peak_group) != 2L) {
  stop("The specified cutoff did not produce both melanoma peak groups.")
}

peak_colors <- c("Left peak" = "#4CAF50", "Right peak" = "#E64B35")
peak_medians <- melanoma |>
  group_by(.data$peak_group) |>
  summarise(median_expression = median(.data$TNFRSF1A), .groups = "drop")

peak_plot <- ggplot(
  melanoma,
  aes(.data$TNFRSF1A, fill = .data$peak_group, color = .data$peak_group)
) +
  geom_histogram(
    aes(y = after_stat(density)), bins = 18, alpha = 0.22,
    position = "identity", color = NA
  ) +
  geom_density(linewidth = 1.5, alpha = 0.85, adjust = 1.2) +
  geom_vline(
    data = peak_medians,
    aes(xintercept = .data$median_expression, color = .data$peak_group),
    linetype = "dashed", linewidth = 1.1, show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = peak_colors) +
  scale_color_manual(values = peak_colors) +
  labs(
    title = "Melanoma: visually defined TNFRSF1A peaks",
    x = paste0(
      "TNFRSF1A expression (",
      gsub("_", " ", liu_units[[1]], fixed = TRUE), ")"
    ),
    y = "Density"
  ) +
  theme_pub()

ggsave(
  file.path(figure_dir, "Supplementary_Figure_S6C_TNFRSF1A_bimodal.png"),
  peak_plot, width = 7.2, height = 5, dpi = 600, bg = "white"
)

purity_data <- melanoma |>
  filter(is.finite(.data$tumor_purity))
if (n_distinct(purity_data$peak_group) != 2L) {
  stop("Tumor purity is not available for both melanoma peak groups.")
}
purity_test <- wilcox.test(tumor_purity ~ peak_group, data = purity_data)
purity_counts <- purity_data |>
  count(.data$peak_group, name = "n")

purity_plot <- ggplot(
  purity_data,
  aes(.data$peak_group, .data$tumor_purity, fill = .data$peak_group)
) +
  geom_boxplot(
    alpha = 0.9, width = 0.58, outlier.shape = NA,
    color = "black", linewidth = 1
  ) +
  geom_jitter(width = 0.08, size = 2.4, alpha = 0.7, color = "black") +
  geom_text(
    data = purity_counts,
    aes(x = .data$peak_group, y = 1.02, label = paste0("n = ", .data$n)),
    inherit.aes = FALSE, size = 4.3
  ) +
  annotate(
    "text", x = 1.5, y = 0.96,
    label = paste0("Wilcoxon p = ", signif(purity_test$p.value, 2)),
    size = 4.6, fontface = "bold"
  ) +
  scale_fill_manual(values = peak_colors) +
  coord_cartesian(ylim = c(0, 1.05), clip = "off") +
  labs(
    title = "Melanoma: tumor purity by TNFRSF1A peak",
    x = NULL, y = "Tumor purity"
  ) +
  theme_pub() +
  theme(legend.position = "none")

ggsave(
  file.path(figure_dir, "Supplementary_Figure_S6C_tumor_purity.png"),
  purity_plot, width = 6, height = 5.2, dpi = 600, bg = "white"
)

peak_summary <- melanoma |>
  group_by(.data$peak_group) |>
  summarise(
    n = n(),
    TNFRSF1A_mean = mean(.data$TNFRSF1A, na.rm = TRUE),
    TNFRSF1A_median = median(.data$TNFRSF1A, na.rm = TRUE),
    purity_mean = mean(.data$tumor_purity, na.rm = TRUE),
    purity_median = median(.data$tumor_purity, na.rm = TRUE),
    IFNG_mean = mean(.data$IFNG, na.rm = TRUE),
    CD8A_mean = mean(.data$CD8A, na.rm = TRUE),
    ICAM1_mean = mean(.data$ICAM1, na.rm = TRUE),
    IRF1_mean = mean(.data$IRF1, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(
  peak_summary,
  file.path(figure_dir, "Supplementary_Figure_S6C_peak_QC_summary.csv")
)
