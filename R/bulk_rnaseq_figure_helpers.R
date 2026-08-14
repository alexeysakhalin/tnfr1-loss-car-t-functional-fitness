# Validation and plotting helpers for the bulk RNA-seq figure adapters.

BULK_RNASEQ_EXPECTED_SYMBOLS <- 46425L

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
      is.finite(adapter$effect) & !is.na(adapter$adjusted_p_value),
    ,
    drop = FALSE
  ]
}

write_bulk_figure_qc <- function(adapter, path, minimum_base_mean = 30) {
  split_adapter <- split(adapter, adapter$condition)
  rows <- lapply(names(split_adapter), function(condition) {
    current <- split_adapter[[condition]]
    tested <- tested_bulk_rows(current, minimum_base_mean)
    triplet <- is.finite(current$effect) & is.finite(current$lfc_se) &
      is.finite(current$wald_statistic)
    expected <- current$effect[triplet] / current$lfc_se[triplet]
    error <- abs(current$wald_statistic[triplet] - expected)
    data.frame(
      condition = condition,
      full_universe_symbols = nrow(current),
      complete_wald_triplets = sum(triplet),
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
