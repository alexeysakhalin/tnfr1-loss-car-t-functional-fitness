# Figure 1B-C: cytokine-induced differential expression in HeLa cells.
#
# Figure 1C can be reproduced from the version-controlled FDR-significant
# tables. Figure 1B requires complete, unfiltered DESeq2 results. The script
# deliberately refuses to draw a conventional volcano plot from a table that
# contains only significant genes.

required_packages <- c(
  "data.table", "R.utils", "dplyr", "ggplot2", "ggrepel", "VennDiagram"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

figure_dir <- file.path("figures", "figure_1")
table_dir <- file.path("results", "figure_1")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

significant_path <- file.path(
  "data", "experimental", "hela_cytokine_significant_differential_expression.tsv.gz"
)
unfiltered_path <- file.path(
  "data", "experimental", "hela_cytokine_unfiltered_differential_expression.tsv.gz"
)

required_columns <- c(
  "condition", "gene_symbol", "base_mean",
  "log2_fold_change_treatment_vs_untreated", "lfc_se",
  "wald_statistic_treatment_vs_untreated", "p_value", "adjusted_p_value"
)
expected_conditions <- c("TNF", "IFNG", "TNF_IFNG")

read_de_table <- function(path, label, allow_deseq_missing = FALSE) {
  if (!file.exists(path)) stop(label, " not found: ", path)
  table <- data.table::fread(path, data.table = FALSE)
  missing_columns <- setdiff(required_columns, names(table))
  if (length(missing_columns) > 0L) {
    stop(label, " is missing columns: ", paste(missing_columns, collapse = ", "))
  }
  table <- table[, required_columns]
  if (!setequal(unique(table$condition), expected_conditions)) {
    stop(label, " must contain TNF, IFNG and TNF_IFNG contrasts.")
  }
  if (anyDuplicated(table[c("condition", "gene_symbol")])) {
    stop(label, " contains duplicate condition/gene rows.")
  }
  numeric_columns <- setdiff(required_columns, c("condition", "gene_symbol"))
  if (any(!vapply(table[numeric_columns], is.numeric, logical(1)))) {
    stop(label, " contains non-numeric analysis columns.")
  }
  if (any(is.na(table$base_mean)) || any(!is.finite(table$base_mean)) ||
      any(table$base_mean < 0)) {
    stop(label, " contains invalid base means.")
  }
  effect_columns <- c(
    "log2_fold_change_treatment_vs_untreated", "lfc_se",
    "wald_statistic_treatment_vs_untreated"
  )
  effect_missing <- is.na(table[effect_columns])
  if (any(rowSums(effect_missing) > 0 & rowSums(effect_missing) < 3L)) {
    stop(label, " contains a partially missing LFC/SE/Wald triplet.")
  }
  effect_available <- !is.na(table$log2_fold_change_treatment_vs_untreated)
  if (any(!is.finite(as.matrix(table[effect_available, effect_columns]))) ||
      any(table$lfc_se[effect_available] <= 0)) {
    stop(label, " contains invalid LFC/SE/Wald values.")
  }
  if (!isTRUE(all.equal(
    table$log2_fold_change_treatment_vs_untreated[effect_available] /
      table$lfc_se[effect_available],
    table$wald_statistic_treatment_vs_untreated[effect_available],
    tolerance = 1e-8,
    check.attributes = FALSE
  ))) {
    stop(label, " contains an inconsistent LFC/SE/Wald relationship.")
  }
  p_available <- !is.na(table$p_value)
  padj_available <- !is.na(table$adjusted_p_value)
  if (any(p_available & !effect_available) || any(padj_available & !p_available)) {
    stop(label, " contains an invalid DESeq2 p-value missingness pattern.")
  }
  if (any(!is.finite(table$p_value[p_available])) ||
      any(table$p_value[p_available] < 0 | table$p_value[p_available] > 1) ||
      any(!is.finite(table$adjusted_p_value[padj_available])) ||
      any(table$adjusted_p_value[padj_available] < 0 |
          table$adjusted_p_value[padj_available] > 1)) {
    stop(label, " contains invalid p-values.")
  }
  if (!allow_deseq_missing && any(is.na(table[numeric_columns]))) {
    stop(label, " must not contain missing analysis values.")
  }
  table$gene_symbol <- toupper(trimws(table$gene_symbol))
  if (any(table$gene_symbol == "")) stop(label, " contains blank gene symbols.")
  table
}

significant <- read_de_table(significant_path, "Significant-gene DE table")
if (any(significant$adjusted_p_value >= 0.05)) {
  stop("The version-controlled significant-gene table must contain only FDR < 0.05 rows.")
}

eligible_significant <- significant |>
  filter(base_mean >= 30, log2_fold_change_treatment_vs_untreated > 1)

up_sets <- lapply(expected_conditions, function(condition_name) {
  sort(unique(eligible_significant$gene_symbol[
    eligible_significant$condition == condition_name
  ]))
})
names(up_sets) <- expected_conditions

venn_counts <- data.frame(
  condition = c("TNF", "IFNG", "TNF_IFNG", "TNF_IFNG_common"),
  n = c(
    length(up_sets$TNF),
    length(up_sets$IFNG),
    length(up_sets$TNF_IFNG),
    length(Reduce(intersect, up_sets))
  )
)
data.table::fwrite(
  venn_counts,
  file.path(table_dir, "Figure_1C_upregulated_gene_counts.tsv"),
  sep = "\t"
)

common_upregulated <- sort(Reduce(intersect, up_sets))
common_table <- eligible_significant |>
  filter(gene_symbol %in% common_upregulated) |>
  select(
    condition, gene_symbol, base_mean,
    log2_fold_change_treatment_vs_untreated, adjusted_p_value
  ) |>
  arrange(gene_symbol, factor(condition, levels = expected_conditions))
data.table::fwrite(
  common_table,
  file.path(table_dir, "Figure_1C_common_upregulated_genes.tsv"),
  sep = "\t"
)

n12 <- length(intersect(up_sets$TNF, up_sets$IFNG))
n13 <- length(intersect(up_sets$TNF, up_sets$TNF_IFNG))
n23 <- length(intersect(up_sets$IFNG, up_sets$TNF_IFNG))
n123 <- length(common_upregulated)

venn_grob <- VennDiagram::draw.triple.venn(
  area1 = length(up_sets$TNF),
  area2 = length(up_sets$IFNG),
  area3 = length(up_sets$TNF_IFNG),
  n12 = n12,
  n13 = n13,
  n23 = n23,
  n123 = n123,
  category = c("TNF", "IFNγ", "TNF + IFNγ"),
  fill = c("#6FC7CF", "#1CC5FE", "#FBA27D"),
  alpha = 0.75,
  cex = 1.5,
  cat.cex = 1.4,
  lwd = 1.5,
  scaled = FALSE,
  ind = FALSE
)

save_venn <- function(path, device) {
  if (identical(device, "png")) {
    grDevices::png(path, width = 6, height = 6, units = "in", res = 600)
  } else {
    grDevices::tiff(
      path, width = 6, height = 6, units = "in", res = 600,
      compression = "lzw"
    )
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.draw(venn_grob)
}
save_venn(file.path(figure_dir, "Figure_1C_upregulated_overlap.png"), "png")
save_venn(file.path(figure_dir, "Figure_1C_upregulated_overlap.tiff"), "tiff")
writeLines(capture.output(sessionInfo()), file.path(table_dir, "sessionInfo.txt"))

if (!file.exists(unfiltered_path)) {
  message(
    "Figure 1C was generated. Figure 1B remains locked because complete ",
    "unfiltered DESeq2 results are absent at ", unfiltered_path, "."
  )
} else {

unfiltered <- read_de_table(
  unfiltered_path, "Unfiltered DE table", allow_deseq_missing = TRUE
)
unfiltered_gene_sets <- lapply(expected_conditions, function(condition_name) {
  sort(unfiltered$gene_symbol[unfiltered$condition == condition_name])
})
if (!all(vapply(
  unfiltered_gene_sets[-1], identical, logical(1), unfiltered_gene_sets[[1]]
))) {
  stop("The unfiltered contrasts do not share one complete gene universe.")
}
if (length(unfiltered_gene_sets[[1]]) < 15000L) {
  stop("The unfiltered gene universe is unexpectedly small (<15,000 genes).")
}

comparison_columns <- c(
  "base_mean", "log2_fold_change_treatment_vs_untreated", "lfc_se",
  "wald_statistic_treatment_vs_untreated", "p_value", "adjusted_p_value"
)
significant_check <- significant |>
  inner_join(
    unfiltered,
    by = c("condition", "gene_symbol"),
    suffix = c("_significant", "_unfiltered")
  )
if (nrow(significant_check) != nrow(significant)) {
  stop("The unfiltered table is missing one or more canonical significant-gene rows.")
}
for (column_name in comparison_columns) {
  significant_values <- significant_check[[paste0(column_name, "_significant")]]
  unfiltered_values <- significant_check[[paste0(column_name, "_unfiltered")]]
  if (!isTRUE(all.equal(
    significant_values, unfiltered_values,
    tolerance = 1e-10, check.attributes = FALSE
  ))) {
    stop(
      "The unfiltered and canonical significant tables disagree in ", column_name, "."
    )
  }
}
condition_qc <- unfiltered |>
  group_by(condition) |>
  summarise(
    n_genes = n(),
    n_adjusted_p_available = sum(!is.na(adjusted_p_value)),
    n_not_significant = sum(adjusted_p_value >= 0.05, na.rm = TRUE),
    .groups = "drop"
  )
if (any(condition_qc$n_not_significant < 1000L)) {
  stop("The purported unfiltered table has fewer than 1,000 non-significant genes in a contrast.")
}
data.table::fwrite(
  condition_qc,
  file.path(table_dir, "Figure_1B_input_qc.tsv"),
  sep = "\t"
)

label_genes <- c(
  "CASP3", "CASP7", "CASP8", "CASP9", "BAX", "BAK1", "BCL2", "FAS",
  "APAF1", "GSDMD", "GSDME", "CASP1", "CASP4", "CASP5", "AIM2",
  "NLRP3", "RIPK1", "RIPK3", "MLKL", "ICAM1", "IRF1"
)
condition_titles <- c(TNF = "TNF", IFNG = "IFNγ", TNF_IFNG = "TNF + IFNγ")
condition_filenames <- c(TNF = "TNF", IFNG = "IFNg", TNF_IFNG = "TNF_IFNg")

make_volcano <- function(condition_name) {
  plot_data <- unfiltered |>
    filter(
      condition == condition_name,
      base_mean >= 30,
      !is.na(log2_fold_change_treatment_vs_untreated),
      !is.na(adjusted_p_value)
    ) |>
    mutate(
      adjusted_p_plot = pmax(adjusted_p_value, .Machine$double.xmin),
      significance = case_when(
        log2_fold_change_treatment_vs_untreated > 1 & adjusted_p_value < 0.05 ~
          "Upregulated",
        log2_fold_change_treatment_vs_untreated < -1 & adjusted_p_value < 0.05 ~
          "Downregulated",
        TRUE ~ "Not significant"
      )
    )
  labels <- plot_data |>
    filter(gene_symbol %in% label_genes)
  plot <- ggplot(
    plot_data,
    aes(
      x = log2_fold_change_treatment_vs_untreated,
      y = -log10(adjusted_p_plot)
    )
  ) +
    geom_point(aes(color = significance), size = 1.1, alpha = 0.75) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.35) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.35) +
    ggrepel::geom_text_repel(
      data = labels,
      aes(label = gene_symbol),
      seed = 42,
      size = 3.2,
      min.segment.length = 0,
      max.overlaps = Inf
    ) +
    scale_color_manual(
      values = c(
        "Upregulated" = "#D1495B",
        "Downregulated" = "#2B6CB0",
        "Not significant" = "grey75"
      ),
      breaks = c("Upregulated", "Downregulated", "Not significant")
    ) +
    labs(
      title = condition_titles[[condition_name]],
      x = expression(log[2] * " fold change (treatment vs untreated)"),
      y = expression(-log[10] * " adjusted p-value"),
      color = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5))

  base_name <- paste0(
    "Figure_1B_", condition_filenames[[condition_name]], "_volcano"
  )
  ggsave(
    file.path(figure_dir, paste0(base_name, ".png")),
    plot, width = 7, height = 5.5, units = "in", dpi = 600, bg = "white"
  )
  ggsave(
    file.path(figure_dir, paste0(base_name, ".tiff")),
    plot, width = 7, height = 5.5, units = "in", dpi = 600,
    compression = "lzw", bg = "white"
  )
}

invisible(lapply(expected_conditions, make_volcano))
message("Figure 1B-C outputs written to ", figure_dir)
}
