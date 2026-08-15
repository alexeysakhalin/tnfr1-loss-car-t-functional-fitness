# Figure 2B and Supplementary Figure S2D: TNFR1-KO1 versus WT.
#
# The canonical source calls TNFR1-KO1 "T6". The adapter used here is generated
# from the integer count matrix by scripts/run_bulk_rnaseq_pydeseq2.py. Positive
# log2 fold change means higher expression in TNFR1-KO1 than in WT within the
# indicated treatment stratum; these are not genotype-by-treatment interactions.

required_packages <- c("ggplot2", "ggrepel", "VennDiagram")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "), ".")
}

source(file.path("R", "bulk_rnaseq_figure_helpers.R"))

input_path <- file.path(
  "data", "experimental", "bulk_rnaseq", "derived",
  "figure_2b_s2d_tnfr1_ko1_vs_wt_matched_treatments.unfiltered.tsv.gz"
)
expected_conditions <- c("control", "IFNG", "TNF", "TNF_IFNG")
effect_column <- "log2_fold_change_ko1_vs_wt"
statistic_column <- "wald_statistic_ko1_vs_wt"

ko1_vs_wt <- read_bulk_figure_adapter(
  path = input_path,
  expected_conditions = expected_conditions,
  effect_column = effect_column,
  statistic_column = statistic_column
)

figure_dir <- file.path("figures", "figure_2")
result_dir <- file.path("results", "figure_2")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
write_bulk_figure_qc(
  ko1_vs_wt,
  file.path(result_dir, "Figure_2B_S2D_input_qc.tsv")
)

condition_display <- c(
  control = "Untreated",
  IFNG = "IFNγ",
  TNF = "TNF",
  TNF_IFNG = "TNF + IFNγ"
)
if (!identical(names(condition_display), expected_conditions)) {
  stop("Internal condition-display mapping is incomplete.")
}

figure_2b_label_genes <- c(
  "ICAM1", "MLKL", "GSDME", "IRF1"
)
publication_condition_order <- c("TNF_IFNG", "IFNG", "TNF")
ko1_prepared <- setNames(
  lapply(publication_condition_order, function(condition) {
    condition_rows <- ko1_vs_wt[
      ko1_vs_wt$condition == condition,
      ,
      drop = FALSE
    ]
    prepare_bulk_volcano_rows(condition_rows, figure_2b_label_genes)
  }),
  publication_condition_order
)
ko1_common_x_limits <- symmetric_bulk_volcano_x_limits(ko1_prepared)
write_bulk_volcano_output_contract(
  prepared_by_condition = ko1_prepared,
  condition_order = publication_condition_order,
  label_genes = figure_2b_label_genes,
  common_x_limits = ko1_common_x_limits,
  path = file.path(result_dir, "Figure_2B_volcano_output_contract.tsv")
)

make_ko1_volcano_plot <- function(condition, panel_title, compact = FALSE) {
  prepared <- ko1_prepared[[condition]]
  if (is.null(prepared)) {
    stop("No prepared volcano data for condition ", condition, ".")
  }
  genes <- prepared$genes
  label_data <- prepared$labels
  base_size <- if (compact) 8.5 else 14
  label_size <- if (compact) BULK_VOLCANO_COMPACT_LABEL_SIZE_MM else 4
  point_size <- if (compact) 0.70 else 1.8
  x_label <- if (compact) NULL else {
    "Log2 fold change (TNFR1-KO1 vs WT)"
  }
  y_label <- if (compact) NULL else {
    "-Log10 adjusted P (capped at 300)"
  }

  plot <- ggplot2::ggplot(
    genes,
    ggplot2::aes(x = effect, y = minus_log10_adjusted_p)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = status),
      size = point_size,
      alpha = 0.80
    ) +
    ggplot2::scale_color_manual(
      values = BULK_VOLCANO_STATUS_COLORS,
      drop = FALSE,
      name = "Status"
    ) +
    ggplot2::geom_vline(
      xintercept = c(-1, 1), linetype = "dashed", color = "black"
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(0.05), linetype = "dashed", color = "black"
    ) +
    ggplot2::scale_y_continuous(
      breaks = c(0, 50, 100, 150, 200, 250, 300),
      expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    ggplot2::coord_cartesian(
      xlim = ko1_common_x_limits,
      ylim = c(0, BULK_VOLCANO_Y_CAP + 5),
      clip = "off"
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::labs(
      title = panel_title,
      x = x_label,
      y = y_label
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = if (compact) 9 else ggplot2::rel(1.2),
        face = if (compact) "bold" else "plain"
      ),
      axis.text = ggplot2::element_text(size = if (compact) 8 else ggplot2::rel(0.8)),
      axis.title = ggplot2::element_text(size = if (compact) 8.5 else ggplot2::rel(1)),
      legend.position = if (compact) "none" else "right",
      plot.margin = if (compact) {
        grid::unit(c(4, 3, 2, 2), "pt")
      } else {
        ggplot2::margin()
      }
    ) +
    ggrepel::geom_label_repel(
      data = label_data,
      ggplot2::aes(label = gene_symbol),
      color = "black",
      size = label_size,
      box.padding = if (compact) 0.15 else 0.35,
      point.padding = if (compact) 0.10 else 0.25,
      label.padding = grid::unit(if (compact) 0.08 else 0.15, "lines"),
      label.size = 0,
      fill = scales::alpha("white", 0.9),
      segment.color = "black",
      segment.size = if (compact) 0.18 else 0.3,
      min.segment.length = 0,
      max.overlaps = Inf,
      max.time = 2,
      seed = 42
    )
  plot
}

ko1_individual_plots <- list(
  TNF_IFNG = make_ko1_volcano_plot(
    "TNF_IFNG",
    expression("TNF + IFN" * gamma)
  ),
  IFNG = make_ko1_volcano_plot("IFNG", expression("IFN" * gamma)),
  TNF = make_ko1_volcano_plot("TNF", "TNF")
)
ko1_individual_stems <- c(
  TNF_IFNG = "Figure_2B_TNF_IFNg",
  IFNG = "Figure_2B_IFNg",
  TNF = "Figure_2B_TNF"
)
for (condition in publication_condition_order) {
  save_figure_pair(
    ko1_individual_plots[[condition]],
    file.path(figure_dir, ko1_individual_stems[[condition]]),
    width = 8,
    height = 6
  )
}

ko1_triptych_plots <- list(
  make_ko1_volcano_plot("TNF_IFNG", expression("TNF + IFN" * gamma), compact = TRUE),
  make_ko1_volcano_plot("IFNG", expression("IFN" * gamma), compact = TRUE),
  make_ko1_volcano_plot("TNF", "TNF", compact = TRUE)
)
save_bulk_volcano_triptych(
  plots = ko1_triptych_plots,
  stem = file.path(figure_dir, "Figure_2B_triptych"),
  shared_x_label = "Log2 fold change (TNFR1-KO1 vs WT)",
  shared_y_label = "-Log10 adjusted P (capped at 300)"
)

# Supplementary Figure S2D is restricted to the three cytokine-treated strata.
figure_conditions <- c("TNF", "IFNG", "TNF_IFNG")
tested_by_condition <- setNames(
  lapply(figure_conditions, function(condition) {
    tested_bulk_rows(ko1_vs_wt[ko1_vs_wt$condition == condition, , drop = FALSE])
  }),
  figure_conditions
)
downregulated_symbols <- lapply(tested_by_condition, function(frame) {
  unique(frame$gene_symbol[
    frame$effect < -1 & frame$adjusted_p_value < 0.05
  ])
})

down_tnf <- downregulated_symbols[["TNF"]]
down_ifng <- downregulated_symbols[["IFNG"]]
down_tnf_ifng <- downregulated_symbols[["TNF_IFNG"]]

venn_plot <- VennDiagram::draw.triple.venn(
  area1 = length(down_tnf),
  area2 = length(down_ifng),
  area3 = length(down_tnf_ifng),
  n12 = length(intersect(down_tnf, down_ifng)),
  n23 = length(intersect(down_ifng, down_tnf_ifng)),
  n13 = length(intersect(down_tnf, down_tnf_ifng)),
  n123 = length(Reduce(intersect, list(down_tnf, down_ifng, down_tnf_ifng))),
  category = c("TNF", "IFNγ", "TNF + IFNγ"),
  fill = c("#6FC7CF", "#1CC5FE", "#FBA27D"),
  alpha = 0.75,
  cex = 2,
  fontface = "bold",
  cat.cex = 2,
  cat.fontface = "bold",
  lwd = 2,
  scaled = FALSE,
  ind = FALSE
)

save_venn_pair <- function(venn, stem) {
  grDevices::tiff(
    paste0(stem, ".tiff"),
    width = 6,
    height = 6,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  grid::grid.newpage()
  grid::grid.draw(venn)
  grDevices::dev.off()

  grDevices::png(
    paste0(stem, ".png"),
    width = 6,
    height = 6,
    units = "in",
    res = 600,
    bg = "white"
  )
  grid::grid.newpage()
  grid::grid.draw(venn)
  grDevices::dev.off()
}
save_venn_pair(
  venn_plot,
  file.path(figure_dir, "Supplementary_Figure_S2D_downregulated_overlap")
)

common_downregulated <- sort(
  Reduce(intersect, list(down_tnf, down_ifng, down_tnf_ifng))
)
common_stats <- data.frame(
  gene_symbol = common_downregulated,
  stringsAsFactors = FALSE
)
for (condition in figure_conditions) {
  current <- tested_by_condition[[condition]]
  current <- current[current$gene_symbol %in% common_downregulated, , drop = FALSE]
  suffix <- switch(
    condition,
    TNF = "TNF",
    IFNG = "IFNg",
    TNF_IFNG = "TNF_IFNg"
  )
  condition_stats <- current[
    , c("gene_symbol", "effect", "adjusted_p_value", "base_mean"), drop = FALSE
  ]
  names(condition_stats)[-1] <- c(
    paste0("log2_fold_change_KO1_vs_WT_", suffix),
    paste0("adjusted_p_value_", suffix),
    paste0("base_mean_", suffix)
  )
  common_stats <- merge(
    common_stats,
    condition_stats,
    by = "gene_symbol",
    all.x = TRUE,
    sort = FALSE
  )
  common_stats <- common_stats[
    match(common_downregulated, common_stats$gene_symbol),
    ,
    drop = FALSE
  ]
}

utils::write.table(
  data.frame(gene_symbol = common_downregulated, stringsAsFactors = FALSE),
  file = file.path(
    result_dir,
    "Supplementary_Figure_S2D_common_downregulated_genes.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "NA"
)
utils::write.table(
  common_stats,
  file = file.path(
    result_dir,
    "Supplementary_Figure_S2D_common_downregulated_gene_statistics.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "NA"
)
