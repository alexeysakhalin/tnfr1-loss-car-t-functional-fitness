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

label_genes <- c(
  "CASP3", "CASP7", "CASP8", "CASP9", "BAX", "BAK1", "BCL2", "FAS",
  "APAF1", "GSDMD", "GSDME", "CASP1", "CASP4", "CASP5", "AIM2",
  "NLRP3", "RIPK1", "RIPK3", "MLKL", "ICAM1", "IRF1"
)

make_ko1_volcano <- function(condition, panel_title, output_stem) {
  condition_rows <- ko1_vs_wt[ko1_vs_wt$condition == condition, , drop = FALSE]
  genes <- tested_bulk_rows(condition_rows)
  if (!nrow(genes)) {
    stop("No tested genes remain for the ", condition, " volcano plot.")
  }

  genes$adjusted_p_value_plot <- pmax(genes$adjusted_p_value, 1e-300)
  genes$minus_log10_adjusted_p <- -log10(genes$adjusted_p_value_plot)
  genes$status <- factor(
    ifelse(
      genes$effect > 1 & genes$adjusted_p_value < 0.05,
      "Upregulated",
      ifelse(
        genes$effect < -1 & genes$adjusted_p_value < 0.05,
        "Downregulated",
        "Not significant"
      )
    ),
    levels = c("Upregulated", "Downregulated", "Not significant")
  )

  label_data <- genes[genes$gene_symbol_key %in% label_genes, , drop = FALSE]
  y_limit <- max(20, ceiling(max(genes$minus_log10_adjusted_p, na.rm = TRUE)) + 2)

  plot <- ggplot2::ggplot(
    genes,
    ggplot2::aes(x = effect, y = minus_log10_adjusted_p)
  ) +
    ggplot2::geom_point(ggplot2::aes(color = status), size = 1.8, alpha = 0.80) +
    ggplot2::scale_color_manual(
      values = c(
        "Upregulated" = "#E64B35",
        "Downregulated" = "#3182BD",
        "Not significant" = "grey70"
      ),
      drop = FALSE,
      name = "Status"
    ) +
    ggplot2::geom_vline(
      xintercept = c(-1, 1), linetype = "dashed", color = "black"
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(0.05), linetype = "dashed", color = "black"
    ) +
    ggplot2::coord_cartesian(xlim = c(-15, 15), ylim = c(0, y_limit), clip = "on") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      title = panel_title,
      x = "Log2 fold change (TNFR1-KO1 vs WT)",
      y = "-Log10 adjusted p-value"
    ) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggrepel::geom_label_repel(
      data = label_data,
      ggplot2::aes(label = gene_symbol),
      color = "black",
      size = 4,
      box.padding = 0.35,
      point.padding = 0.25,
      label.padding = grid::unit(0.15, "lines"),
      label.size = 0,
      fill = scales::alpha("white", 0.9),
      segment.color = "black",
      segment.size = 0.3,
      min.segment.length = 0,
      max.overlaps = Inf,
      seed = 42
    )

  save_figure_pair(
    plot,
    file.path(figure_dir, output_stem),
    width = 8,
    height = 6
  )
  invisible(plot)
}

make_ko1_volcano(
  "TNF_IFNG",
  expression("TNF + IFN" * gamma),
  "Figure_2B_TNF_IFNg"
)
make_ko1_volcano("TNF", "TNF", "Figure_2B_TNF")
make_ko1_volcano("IFNG", expression("IFN" * gamma), "Figure_2B_IFNg")

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
