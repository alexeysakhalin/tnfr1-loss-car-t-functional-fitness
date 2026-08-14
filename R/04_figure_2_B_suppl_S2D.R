# Figure 2B and Supplementary Figure S2D: matched T6-versus-WT RNA-seq.
#
# T6 is the internal identifier present in the supplied design workbook. The
# script does not relabel T6 as KO1 or KO2 because that mapping is not encoded
# in the data and must be confirmed from the experimental sample key.

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

figure_dir <- file.path("figures", "figure_2")
table_dir <- file.path("results", "figure_2")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

input_paths <- c(
  UNTREATED = file.path(
    "data", "experimental", "hela_t6_vs_wt_untreated_differential_expression.tsv.gz"
  ),
  IFN = file.path(
    "data", "experimental", "hela_t6_vs_wt_ifn_differential_expression.tsv.gz"
  ),
  TNF = file.path(
    "data", "experimental", "hela_t6_vs_wt_tnf_differential_expression.tsv.gz"
  ),
  TI = file.path(
    "data", "experimental", "hela_t6_vs_wt_ti_differential_expression.tsv.gz"
  )
)

required_columns <- c(
  "condition", "gene_symbol", "base_mean", "log2_fold_change_t6_vs_wt",
  "lfc_se", "wald_statistic_t6_vs_wt", "p_value", "adjusted_p_value"
)

read_condition <- function(condition_name, path) {
  if (!file.exists(path)) stop("Missing matched DE table: ", path)
  table <- data.table::fread(path, data.table = FALSE)
  missing_columns <- setdiff(required_columns, names(table))
  if (length(missing_columns) > 0L) {
    stop(path, " is missing columns: ", paste(missing_columns, collapse = ", "))
  }
  table <- table[, required_columns]
  if (!identical(unique(table$condition), condition_name)) {
    stop("Condition label mismatch in ", path)
  }
  table$gene_symbol <- toupper(trimws(table$gene_symbol))
  if (any(table$gene_symbol == "") || anyDuplicated(table$gene_symbol)) {
    stop("Blank or duplicate gene symbols in ", path)
  }
  if (!is.numeric(table$base_mean) || any(!is.finite(table$base_mean)) ||
      any(table$base_mean < 0)) {
    stop("Invalid base_mean values in ", path)
  }
  effect_columns <- c(
    "log2_fold_change_t6_vs_wt", "lfc_se", "wald_statistic_t6_vs_wt"
  )
  effect_missing <- is.na(table[effect_columns])
  if (any(rowSums(effect_missing) > 0 & rowSums(effect_missing) < 3L)) {
    stop("Partially missing LFC/SE/Wald triplet in ", path)
  }
  effect_available <- !is.na(table$log2_fold_change_t6_vs_wt)
  if (any(!is.finite(as.matrix(table[effect_available, effect_columns])))) {
    stop("Non-finite LFC/SE/Wald values in ", path)
  }
  if (any(table$lfc_se[effect_available] <= 0)) {
    stop("Non-positive LFC standard error in ", path)
  }
  if (!isTRUE(all.equal(
      table$log2_fold_change_t6_vs_wt[effect_available] /
        table$lfc_se[effect_available],
      table$wald_statistic_t6_vs_wt[effect_available],
      tolerance = 1e-8,
      check.attributes = FALSE
    ))) {
    stop("LFC/SE/Wald relationship is inconsistent in ", path)
  }

  # DESeq2 can report a finite effect and p-value with padj=NA after independent
  # filtering, or leave p and padj missing for rows that were not tested. These
  # are valid missingness patterns and must not be coerced to significance.
  p_available <- !is.na(table$p_value)
  padj_available <- !is.na(table$adjusted_p_value)
  if (any(p_available & !effect_available) || any(padj_available & !p_available)) {
    stop("Invalid p-value missingness pattern in ", path)
  }
  if (any(!is.finite(table$p_value[p_available])) ||
      any(table$p_value[p_available] < 0 | table$p_value[p_available] > 1) ||
      any(!is.finite(table$adjusted_p_value[padj_available])) ||
      any(table$adjusted_p_value[padj_available] < 0 |
          table$adjusted_p_value[padj_available] > 1)) {
    stop("Invalid DE p-value in ", path)
  }
  table
}

tables <- Map(read_condition, names(input_paths), unname(input_paths))
names(tables) <- names(input_paths)
gene_sets <- lapply(tables, function(table) sort(table$gene_symbol))
if (!all(vapply(gene_sets[-1], identical, logical(1), gene_sets[[1]]))) {
  stop("The four matched contrasts do not have the same gene universe.")
}

input_qc <- data.frame(
  condition = names(tables),
  n_genes = vapply(tables, nrow, integer(1)),
  n_effect_estimates = vapply(
    tables,
    function(table) sum(!is.na(table$log2_fold_change_t6_vs_wt)),
    integer(1)
  ),
  n_p_values = vapply(
    tables, function(table) sum(!is.na(table$p_value)), integer(1)
  ),
  n_adjusted_p_values = vapply(
    tables, function(table) sum(!is.na(table$adjusted_p_value)), integer(1)
  )
)
data.table::fwrite(input_qc, file.path(table_dir, "matched_DE_input_qc.tsv"), sep = "\t")

label_genes <- c(
  "CASP3", "CASP7", "CASP8", "CASP9", "BAX", "BAK1", "BCL2", "FAS",
  "APAF1", "GSDMD", "GSDME", "CASP1", "CASP4", "CASP5", "AIM2",
  "NLRP3", "RIPK1", "RIPK3", "MLKL", "ICAM1", "IRF1"
)
condition_titles <- c(TNF = "TNF", IFN = "IFNγ", TI = "TNF + IFNγ")
condition_filenames <- c(TNF = "TNF", IFN = "IFNg", TI = "TNF_IFNg")

make_volcano <- function(condition_name) {
  plot_data <- tables[[condition_name]] |>
    filter(
      base_mean >= 30,
      !is.na(log2_fold_change_t6_vs_wt),
      !is.na(adjusted_p_value)
    ) |>
    mutate(
      adjusted_p_plot = pmax(adjusted_p_value, .Machine$double.xmin),
      significance = case_when(
        log2_fold_change_t6_vs_wt > 1 & adjusted_p_value < 0.05 ~ "Upregulated",
        log2_fold_change_t6_vs_wt < -1 & adjusted_p_value < 0.05 ~ "Downregulated",
        TRUE ~ "Not significant"
      )
    )
  labels <- plot_data |>
    filter(gene_symbol %in% label_genes)

  plot <- ggplot(
    plot_data,
    aes(x = log2_fold_change_t6_vs_wt, y = -log10(adjusted_p_plot))
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
      subtitle = "T6 versus WT within treatment condition",
      x = expression(log[2] * " fold change (T6 vs WT)"),
      y = expression(-log[10] * " adjusted p-value"),
      color = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )

  output_base <- paste0(
    "Figure_2B_", condition_filenames[[condition_name]], "_T6_vs_WT_volcano"
  )
  ggsave(
    file.path(figure_dir, paste0(output_base, ".png")),
    plot, width = 7, height = 5.5, units = "in", dpi = 600, bg = "white"
  )
  ggsave(
    file.path(figure_dir, paste0(output_base, ".tiff")),
    plot, width = 7, height = 5.5, units = "in", dpi = 600,
    compression = "lzw", bg = "white"
  )
}
invisible(lapply(c("TNF", "IFN", "TI"), make_volcano))

down_sets <- lapply(c("TNF", "IFN", "TI"), function(condition_name) {
  table <- tables[[condition_name]]
  sort(unique(table$gene_symbol[
    table$base_mean >= 30 &
      !is.na(table$log2_fold_change_t6_vs_wt) &
      !is.na(table$adjusted_p_value) &
      table$log2_fold_change_t6_vs_wt < -1 &
      table$adjusted_p_value < 0.05
  ]))
})
names(down_sets) <- c("TNF", "IFN", "TI")

n12 <- length(intersect(down_sets$TNF, down_sets$IFN))
n13 <- length(intersect(down_sets$TNF, down_sets$TI))
n23 <- length(intersect(down_sets$IFN, down_sets$TI))
common_down <- sort(Reduce(intersect, down_sets))
n123 <- length(common_down)

venn_grob <- VennDiagram::draw.triple.venn(
  area1 = length(down_sets$TNF),
  area2 = length(down_sets$IFN),
  area3 = length(down_sets$TI),
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
save_venn(
  file.path(figure_dir, "Supplementary_Figure_S2D_downregulated_overlap.png"),
  "png"
)
save_venn(
  file.path(figure_dir, "Supplementary_Figure_S2D_downregulated_overlap.tiff"),
  "tiff"
)

common_table <- bind_rows(lapply(names(down_sets), function(condition_name) {
  tables[[condition_name]] |>
    filter(gene_symbol %in% common_down) |>
    transmute(
      condition = condition_name,
      gene_symbol,
      base_mean,
      log2_fold_change_t6_vs_wt,
      adjusted_p_value
    )
})) |>
  arrange(gene_symbol, factor(condition, levels = c("TNF", "IFN", "TI")))
data.table::fwrite(
  common_table,
  file.path(table_dir, "Supplementary_Figure_S2D_common_downregulated_genes.tsv"),
  sep = "\t"
)
writeLines(capture.output(sessionInfo()), file.path(table_dir, "sessionInfo.txt"))

message("Figure 2B and Supplementary Figure S2D outputs written to ", figure_dir)
