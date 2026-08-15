# Validation and plotting helpers for the bulk RNA-seq figure adapters.

BULK_RNASEQ_EXPECTED_SYMBOLS <- 46425L
BULK_VOLCANO_P_VALUE_FLOOR <- 1e-300
BULK_VOLCANO_Y_CAP <- 300
BULK_VOLCANO_COMPACT_LABEL_SIZE_MM <- 2.82  # 8.0 pt at final 6.5-inch width
BULK_VOLCANO_STATUS_COLORS <- c(
  "Upregulated" = "#E64B35",
  "Downregulated" = "#3182BD",
  "Not significant" = "grey70"
)

read_bulk_figure_adapter <- function(path, expected_conditions,
                                     effect_column, statistic_column) {
  if (!file.exists(path)) {
    stop(
      "Bulk RNA-seq figure adapter not found: ", path, "\n",
      "Run scripts/run_bulk_rnaseq_pydeseq2.py before rendering the figures."
    )
  }

  required_columns <- c(
    "condition", "gene_symbol", "base_mean", effect_column, "lfc_se",
    statistic_column, "p_value", "adjusted_p_value"
  )
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  adapter <- utils::read.delim(
    connection,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    na.strings = "NA",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!identical(names(adapter), required_columns)) {
    stop(
      "Unexpected columns in ", path, ". Expected exactly: ",
      paste(required_columns, collapse = ", "), ". Observed: ",
      paste(names(adapter), collapse = ", "), "."
    )
  }

  adapter$condition <- trimws(as.character(adapter$condition))
  observed_conditions <- sort(unique(adapter$condition))
  if (!identical(observed_conditions, sort(expected_conditions))) {
    stop(
      "Unexpected condition set in ", path, ". Expected ",
      paste(expected_conditions, collapse = ", "), "; observed ",
      paste(observed_conditions, collapse = ", "), "."
    )
  }

  adapter$gene_symbol <- trimws(as.character(adapter$gene_symbol))
  if (any(is.na(adapter$gene_symbol) | adapter$gene_symbol == "")) {
    stop("Missing or blank gene_symbol values in ", path, ".")
  }

  numeric_columns <- c(
    "base_mean", effect_column, "lfc_se", statistic_column,
    "p_value", "adjusted_p_value"
  )
  for (column in numeric_columns) {
    source_value <- adapter[[column]]
    converted <- suppressWarnings(as.numeric(source_value))
    source_missing <- is.na(source_value) | trimws(as.character(source_value)) == ""
    if (any(!source_missing & is.na(converted))) {
      stop("Non-numeric values found in ", column, " in ", path, ".")
    }
    adapter[[column]] <- converted
  }

  expected_rows <- BULK_RNASEQ_EXPECTED_SYMBOLS * length(expected_conditions)
  if (nrow(adapter) != expected_rows) {
    stop(
      "The adapter is not an unfiltered full-universe table: expected ",
      expected_rows, " rows, observed ", nrow(adapter), "."
    )
  }

  group_counts <- table(factor(adapter$condition, levels = expected_conditions))
  if (any(group_counts != BULK_RNASEQ_EXPECTED_SYMBOLS)) {
    stop(
      "Every condition must contain exactly ", BULK_RNASEQ_EXPECTED_SYMBOLS,
      " symbols; observed ",
      paste(names(group_counts), group_counts, sep = "=", collapse = ", "), "."
    )
  }
  if (anyDuplicated(adapter[c("condition", "gene_symbol")])) {
    stop("Duplicate condition/gene_symbol rows found in ", path, ".")
  }

  reference_symbols <- sort(adapter$gene_symbol[adapter$condition == expected_conditions[[1]]])
  for (condition in expected_conditions[-1]) {
    condition_symbols <- sort(adapter$gene_symbol[adapter$condition == condition])
    if (!identical(condition_symbols, reference_symbols)) {
      stop("Gene-symbol universes differ between conditions in ", path, ".")
    }
  }

  if (any(adapter$base_mean < 0, na.rm = TRUE)) {
    stop("base_mean contains negative values in ", path, ".")
  }
  if (any(adapter$lfc_se <= 0, na.rm = TRUE)) {
    stop("lfc_se must be positive wherever it is reported in ", path, ".")
  }
  for (column in c("p_value", "adjusted_p_value")) {
    value <- adapter[[column]]
    if (any(!is.finite(value[!is.na(value)])) ||
        any(value < 0 | value > 1, na.rm = TRUE)) {
      stop(column, " contains values outside [0, 1] in ", path, ".")
    }
  }
  for (column in c("base_mean", effect_column, "lfc_se", statistic_column)) {
    value <- adapter[[column]]
    if (any(!is.finite(value[!is.na(value)]))) {
      stop(column, " contains non-finite values in ", path, ".")
    }
  }

  complete_triplet <- is.finite(adapter[[effect_column]]) &
    is.finite(adapter$lfc_se) & is.finite(adapter[[statistic_column]])
  if (!any(complete_triplet)) {
    stop("No complete log2-fold-change/standard-error/Wald triplets in ", path, ".")
  }
  triplets_by_condition <- tapply(
    complete_triplet,
    factor(adapter$condition, levels = expected_conditions),
    sum
  )
  if (any(triplets_by_condition == 0)) {
    stop("At least one condition has no complete Wald-test results in ", path, ".")
  }
  expected_statistic <- adapter[[effect_column]][complete_triplet] /
    adapter$lfc_se[complete_triplet]
  observed_statistic <- adapter[[statistic_column]][complete_triplet]
  tolerance <- pmax(1e-8, 1e-6 * abs(expected_statistic))
  if (any(abs(observed_statistic - expected_statistic) > tolerance)) {
    stop(
      "Wald-statistic direction or magnitude is inconsistent with log2FC/lfc_se in ",
      path, "."
    )
  }

  adapter$gene_symbol_key <- toupper(adapter$gene_symbol)
  adapter$effect <- adapter[[effect_column]]
  adapter$wald_statistic <- adapter[[statistic_column]]
  adapter
}

tested_bulk_rows <- function(adapter, minimum_base_mean = 30) {
  adapter[
    !is.na(adapter$base_mean) & adapter$base_mean >= minimum_base_mean &
      is.finite(adapter$effect) & is.finite(adapter$adjusted_p_value),
    ,
    drop = FALSE
  ]
}

plottable_bulk_rows <- function(adapter) {
  adapter[
    is.finite(adapter$effect) & is.finite(adapter$adjusted_p_value),
    ,
    drop = FALSE
  ]
}

prepare_bulk_volcano_rows <- function(adapter, label_genes,
                                       p_value_floor = BULK_VOLCANO_P_VALUE_FLOOR,
                                       y_cap = BULK_VOLCANO_Y_CAP) {
  genes <- plottable_bulk_rows(adapter)
  if (!nrow(genes)) {
    stop("No modelled genes with finite adjusted p-values remain for plotting.")
  }
  if (!is.numeric(p_value_floor) || length(p_value_floor) != 1L ||
      !is.finite(p_value_floor) || p_value_floor <= 0 || p_value_floor >= 1) {
    stop("p_value_floor must be one finite number strictly between 0 and 1.")
  }
  if (!is.numeric(y_cap) || length(y_cap) != 1L ||
      !is.finite(y_cap) || y_cap <= 0) {
    stop("y_cap must be one positive finite number.")
  }

  genes$minus_log10_adjusted_p <- pmin(
    -log10(pmax(genes$adjusted_p_value, p_value_floor)),
    y_cap
  )
  eligible_for_deg <- !is.na(genes$base_mean) & genes$base_mean >= 30
  genes$status <- factor(
    ifelse(
      eligible_for_deg & genes$effect > 1 &
        genes$adjusted_p_value < 0.05,
      "Upregulated",
      ifelse(
        eligible_for_deg & genes$effect < -1 &
          genes$adjusted_p_value < 0.05,
        "Downregulated",
        "Not significant"
      )
    ),
    levels = names(BULK_VOLCANO_STATUS_COLORS)
  )
  label_genes <- unique(toupper(as.character(label_genes)))
  labels <- genes[
    !is.na(genes$base_mean) & genes$base_mean >= 30 &
      genes$gene_symbol_key %in% label_genes,
    ,
    drop = FALSE
  ]
  labels$label_order <- match(labels$gene_symbol_key, label_genes)
  labels <- labels[order(labels$label_order), , drop = FALSE]
  list(genes = genes, labels = labels)
}

symmetric_bulk_volcano_x_limits <- function(prepared_by_condition) {
  effects <- unlist(
    lapply(prepared_by_condition, function(prepared) prepared$genes$effect),
    use.names = FALSE
  )
  if (!length(effects) || any(!is.finite(effects))) {
    stop("Cannot derive a common volcano x range from non-finite effects.")
  }
  outward_limit <- max(1, ceiling(max(abs(effects))))
  c(-outward_limit, outward_limit)
}

write_bulk_volcano_output_contract <- function(
    prepared_by_condition, condition_order, label_genes, common_x_limits, path,
    p_value_floor = BULK_VOLCANO_P_VALUE_FLOOR,
    y_cap = BULK_VOLCANO_Y_CAP) {
  if (!identical(names(prepared_by_condition), condition_order)) {
    stop("Prepared volcano conditions are not in publication panel order.")
  }
  if (length(common_x_limits) != 2L ||
      !all(is.finite(common_x_limits)) ||
      common_x_limits[[1]] >= common_x_limits[[2]]) {
    stop("common_x_limits must contain one increasing finite pair.")
  }
  label_genes <- unique(toupper(as.character(label_genes)))
  rows <- lapply(seq_along(condition_order), function(index) {
    condition <- condition_order[[index]]
    prepared <- prepared_by_condition[[condition]]
    genes <- prepared$genes
    positive_p <- genes$adjusted_p_value[genes$adjusted_p_value > 0]
    finite_uncapped <- -log10(positive_p)
    observed_labels <- label_genes[
      label_genes %in% unique(prepared$labels$gene_symbol_key)
    ]
    data.frame(
      panel_order = index,
      condition = condition,
      plottable_modelled_symbols = nrow(genes),
      minimum_log2_fold_change = min(genes$effect),
      maximum_log2_fold_change = max(genes$effect),
      maximum_absolute_log2_fold_change = max(abs(genes$effect)),
      common_x_min = common_x_limits[[1]],
      common_x_max = common_x_limits[[2]],
      adjusted_p_value_floor = p_value_floor,
      displayed_y_cap = y_cap,
      points_at_displayed_y_cap = sum(
        genes$adjusted_p_value <= p_value_floor
      ),
      maximum_finite_uncapped_minus_log10_adjusted_p = if (length(finite_uncapped)) {
        max(finite_uncapped)
      } else {
        NA_real_
      },
      requested_label_count = length(label_genes),
      plotted_label_count = length(observed_labels),
      plotted_label_genes = paste(observed_labels, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  contract <- do.call(rbind, rows)
  serialized_contract <- contract
  double_columns <- names(serialized_contract)[
    vapply(serialized_contract, is.double, logical(1))
  ]
  for (column in double_columns) {
    serialized_contract[[column]] <- ifelse(
      is.na(serialized_contract[[column]]),
      NA_character_,
      sprintf("%.17g", serialized_contract[[column]])
    )
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    serialized_contract,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
  invisible(contract)
}

write_bulk_figure_qc <- function(adapter, path, minimum_base_mean = 30) {
  split_adapter <- split(adapter, adapter$condition)
  rows <- lapply(names(split_adapter), function(condition) {
    current <- split_adapter[[condition]]
    plottable <- plottable_bulk_rows(current)
    tested <- tested_bulk_rows(current, minimum_base_mean)
    triplet <- is.finite(current$effect) & is.finite(current$lfc_se) &
      is.finite(current$wald_statistic)
    expected <- current$effect[triplet] / current$lfc_se[triplet]
    error <- abs(current$wald_statistic[triplet] - expected)
    data.frame(
      condition = condition,
      full_universe_symbols = nrow(current),
      complete_wald_triplets = sum(triplet),
      plottable_modelled_symbols = nrow(plottable),
      tested_base_mean_ge_30 = nrow(tested),
      fdr_below_0_05 = sum(tested$adjusted_p_value < 0.05),
      up_lfc_gt_1_fdr_below_0_05 = sum(
        tested$effect > 1 & tested$adjusted_p_value < 0.05
      ),
      down_lfc_lt_minus_1_fdr_below_0_05 = sum(
        tested$effect < -1 & tested$adjusted_p_value < 0.05
      ),
      maximum_absolute_wald_identity_error = if (length(error)) max(error) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  qc <- do.call(rbind, rows)
  qc <- qc[order(match(qc$condition, unique(adapter$condition))), , drop = FALSE]
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    qc,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
  invisible(qc)
}

save_figure_pair <- function(plot, stem, width, height, dpi = 600) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = paste0(stem, ".tiff"),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    compression = "lzw",
    bg = "white"
  )
  ggplot2::ggsave(
    filename = paste0(stem, ".png"),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
  invisible(plot)
}

draw_bulk_volcano_triptych <- function(plots, shared_x_label, shared_y_label) {
  if (length(plots) != 3L) {
    stop("A publication volcano triptych must contain exactly three panels.")
  }
  grid::grid.newpage()
  layout <- grid::grid.layout(
    nrow = 3,
    ncol = 4,
    widths = grid::unit.c(
      grid::unit(0.25, "in"),
      grid::unit(rep(1, 3), "null")
    ),
    heights = grid::unit.c(
      grid::unit(1, "null"),
      grid::unit(0.23, "in"),
      grid::unit(0.22, "in")
    )
  )
  grid::pushViewport(grid::viewport(layout = layout))
  for (index in seq_along(plots)) {
    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = index + 1L))
    grid::grid.draw(ggplot2::ggplotGrob(plots[[index]]))
    grid::upViewport()
  }

  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid::grid.text(
    shared_y_label,
    rot = 90,
    gp = grid::gpar(fontsize = 8.5)
  )
  grid::upViewport()

  grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 2:4))
  grid::grid.text(shared_x_label, gp = grid::gpar(fontsize = 8.5))
  grid::upViewport()

  legend_x <- c(0.09, 0.40, 0.70)
  grid::pushViewport(grid::viewport(layout.pos.row = 3, layout.pos.col = 2:4))
  grid::grid.points(
    x = grid::unit(legend_x, "npc"),
    y = grid::unit(rep(0.5, 3), "npc"),
    pch = 16,
    size = grid::unit(1.5, "mm"),
    gp = grid::gpar(col = unname(BULK_VOLCANO_STATUS_COLORS))
  )
  grid::grid.text(
    names(BULK_VOLCANO_STATUS_COLORS),
    x = grid::unit(legend_x + 0.022, "npc"),
    y = grid::unit(rep(0.5, 3), "npc"),
    just = "left",
    gp = grid::gpar(fontsize = 8)
  )
  grid::upViewport(2)
  invisible(plots)
}

save_bulk_volcano_triptych <- function(
    plots, stem, shared_x_label, shared_y_label,
    width = 6.5, height = 3, dpi = 600) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  render <- function(open_device) {
    open_device()
    tryCatch(
      draw_bulk_volcano_triptych(plots, shared_x_label, shared_y_label),
      finally = grDevices::dev.off()
    )
  }
  render(function() {
    grDevices::tiff(
      paste0(stem, ".tiff"),
      width = width,
      height = height,
      units = "in",
      res = dpi,
      compression = "lzw",
      bg = "white"
    )
  })
  render(function() {
    grDevices::png(
      paste0(stem, ".png"),
      width = width,
      height = height,
      units = "in",
      res = dpi,
      bg = "white"
    )
  })
  invisible(plots)
}
