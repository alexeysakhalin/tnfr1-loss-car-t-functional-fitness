# Load the analysis-ready cohort tables used by the bulk-transcriptomic figures.
#
# All four patient-level tables are generated locally and excluded from version
# control. The open-cohort preparation script writes data/analysis/open_cohort_*
# and the CheckMate preparation step writes data/processed/checkmate_* from the
# published supplementary workbook.

analysis_expression_columns <- c(
  "cohort_id", "sample_id", "gene_symbol", "expr_value",
  "expression_unit", "rank_percentile"
)

# rank_percentile is defined after mapping source features to gene symbols and
# averaging duplicate features, using all mapped genes in that sample before
# the analysis-gene subset is written.

analysis_metadata_columns <- c(
  "cohort_id", "sample_id", "trial_id", "treatment_arm",
  "os_time_months", "os_event", "pfs_time_months", "pfs_event",
  "tumor_purity"
)

require_table_columns <- function(x, required, object_name) {
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop(
      object_name, " is missing required columns: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  invisible(TRUE)
}

as_numeric_strict <- function(x, column_name, object_name) {
  raw <- trimws(as.character(x))
  missing <- is.na(raw) | raw == ""
  value <- suppressWarnings(as.numeric(raw))
  if (any(!missing & is.na(value))) {
    stop(
      object_name, " column ", column_name,
      " contains a non-numeric value.", call. = FALSE
    )
  }
  value
}

read_analysis_tsv <- function(path, required, object_name) {
  if (!file.exists(path)) {
    stop("Required analysis table is missing: ", path, call. = FALSE)
  }
  out <- readr::read_tsv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    na = c("", "NA", "NaN"),
    show_col_types = FALSE,
    progress = FALSE
  )
  require_table_columns(out, required, object_name)
  out
}

normalise_expression_table <- function(x, object_name) {
  x <- x |>
    dplyr::transmute(
      cohort_id = trimws(as.character(.data$cohort_id)),
      sample_id = trimws(as.character(.data$sample_id)),
      gene_symbol = trimws(as.character(.data$gene_symbol)),
      expr_value = as_numeric_strict(
        .data$expr_value, "expr_value", .env$object_name
      ),
      expression_unit = trimws(as.character(.data$expression_unit)),
      rank_percentile = as_numeric_strict(
        .data$rank_percentile, "rank_percentile", .env$object_name
      )
    )

  invalid_keys <- x |>
    dplyr::filter(
      is.na(.data$cohort_id) | .data$cohort_id == "" |
        is.na(.data$sample_id) | .data$sample_id == "" |
        is.na(.data$gene_symbol) | .data$gene_symbol == ""
    )
  if (nrow(invalid_keys)) {
    stop(object_name, " contains missing cohort, sample or gene identifiers.",
         call. = FALSE)
  }
  if (any(!is.finite(x$expr_value))) {
    stop(object_name, " contains non-finite expression values.", call. = FALSE)
  }
  invalid_rank <- !is.na(x$rank_percentile) &
    (!is.finite(x$rank_percentile) |
       x$rank_percentile < 0 | x$rank_percentile > 1)
  if (any(invalid_rank)) {
    stop(object_name, " contains rank_percentile values outside [0, 1].",
         call. = FALSE)
  }
  if (any(is.na(x$expression_unit) | x$expression_unit == "")) {
    stop(object_name, " contains missing expression units.", call. = FALSE)
  }
  unit_counts <- x |>
    dplyr::distinct(.data$cohort_id, .data$expression_unit) |>
    dplyr::count(.data$cohort_id, name = "n_units")
  if (any(unit_counts$n_units != 1L)) {
    stop(object_name, " must use one declared expression unit per cohort.",
         call. = FALSE)
  }

  duplicate_keys <- x |>
    dplyr::count(.data$cohort_id, .data$sample_id, .data$gene_symbol) |>
    dplyr::filter(.data$n > 1L)
  if (nrow(duplicate_keys)) {
    stop(object_name, " contains duplicate cohort/sample/gene rows.",
         call. = FALSE)
  }
  x
}

normalise_metadata_table <- function(x, object_name) {
  x <- x |>
    dplyr::transmute(
      cohort_id = trimws(as.character(.data$cohort_id)),
      sample_id = trimws(as.character(.data$sample_id)),
      trial_id = trimws(as.character(.data$trial_id)),
      treatment_arm = trimws(as.character(.data$treatment_arm)),
      os_time_months = as_numeric_strict(
        .data$os_time_months, "os_time_months", .env$object_name
      ),
      os_event = as_numeric_strict(
        .data$os_event, "os_event", .env$object_name
      ),
      pfs_time_months = as_numeric_strict(
        .data$pfs_time_months, "pfs_time_months", .env$object_name
      ),
      pfs_event = as_numeric_strict(
        .data$pfs_event, "pfs_event", .env$object_name
      ),
      tumor_purity = as_numeric_strict(
        .data$tumor_purity, "tumor_purity", .env$object_name
      )
    ) |>
    dplyr::mutate(
      trial_id = dplyr::na_if(.data$trial_id, ""),
      treatment_arm = dplyr::na_if(.data$treatment_arm, "")
    )

  invalid_keys <- x |>
    dplyr::filter(
      is.na(.data$cohort_id) | .data$cohort_id == "" |
        is.na(.data$sample_id) | .data$sample_id == ""
    )
  if (nrow(invalid_keys)) {
    stop(object_name, " contains missing cohort or sample identifiers.",
         call. = FALSE)
  }

  duplicate_keys <- x |>
    dplyr::count(.data$cohort_id, .data$sample_id) |>
    dplyr::filter(.data$n > 1L)
  if (nrow(duplicate_keys)) {
    stop(object_name, " contains duplicate cohort/sample rows.", call. = FALSE)
  }

  for (event_col in c("os_event", "pfs_event")) {
    observed <- sort(unique(stats::na.omit(x[[event_col]])))
    if (length(setdiff(observed, c(0, 1)))) {
      stop(object_name, " column ", event_col,
           " must use 1=event and 0=censored.", call. = FALSE)
    }
  }
  for (time_col in c("os_time_months", "pfs_time_months")) {
    if (any(x[[time_col]] < 0, na.rm = TRUE)) {
      stop(object_name, " column ", time_col,
           " contains a negative follow-up time.", call. = FALSE)
    }
  }
  if (any(x$tumor_purity < 0 | x$tumor_purity > 1, na.rm = TRUE)) {
    stop(object_name, " contains tumor_purity values outside [0, 1].",
         call. = FALSE)
  }
  x
}

load_analysis_tables <- function(
    open_expression = file.path(
      "data", "analysis", "open_cohort_selected_expression.tsv.gz"
    ),
    open_metadata = file.path(
      "data", "analysis", "open_cohort_sample_metadata.tsv.gz"
    ),
    checkmate_expression = file.path(
      "data", "processed", "checkmate_selected_expression.tsv.gz"
    ),
    checkmate_metadata = file.path(
      "data", "processed", "checkmate_sample_metadata.tsv.gz"
    )) {
  expression_open <- read_analysis_tsv(
    open_expression, analysis_expression_columns, "open-cohort expression table"
  ) |>
    normalise_expression_table("open-cohort expression table")
  metadata_open <- read_analysis_tsv(
    open_metadata, analysis_metadata_columns, "open-cohort metadata table"
  ) |>
    normalise_metadata_table("open-cohort metadata table")
  expression_checkmate <- read_analysis_tsv(
    checkmate_expression, analysis_expression_columns, "CheckMate expression table"
  ) |>
    normalise_expression_table("CheckMate expression table")
  metadata_checkmate <- read_analysis_tsv(
    checkmate_metadata, analysis_metadata_columns, "CheckMate metadata table"
  ) |>
    normalise_metadata_table("CheckMate metadata table")

  expected_open_cohorts <- c(
    "IMvigor210_BLCA", "SU2C_MARK_NSCLC", "LIU2019_MELANOMA"
  )
  if (!setequal(unique(expression_open$cohort_id), expected_open_cohorts) ||
      !setequal(unique(metadata_open$cohort_id), expected_open_cohorts)) {
    stop(
      "The open-cohort tables must contain exactly: ",
      paste(expected_open_cohorts, collapse = ", "), call. = FALSE
    )
  }
  if (!identical(unique(expression_checkmate$cohort_id), "CHECKMATE_CCRCC") ||
      !identical(unique(metadata_checkmate$cohort_id), "CHECKMATE_CCRCC")) {
    stop("The CheckMate tables must contain only cohort_id=CHECKMATE_CCRCC.",
         call. = FALSE)
  }

  expression <- dplyr::bind_rows(expression_open, expression_checkmate)
  metadata <- dplyr::bind_rows(metadata_open, metadata_checkmate)

  unmatched_expression <- expression |>
    dplyr::anti_join(metadata, by = c("cohort_id", "sample_id"))
  if (nrow(unmatched_expression)) {
    stop("Expression rows without matching sample metadata were found.",
         call. = FALSE)
  }
  unmatched_metadata <- metadata |>
    dplyr::anti_join(
      dplyr::distinct(expression, .data$cohort_id, .data$sample_id),
      by = c("cohort_id", "sample_id")
    )
  if (nrow(unmatched_metadata)) {
    stop("Metadata rows without selected expression were found.",
         call. = FALSE)
  }

  analysis_long <- expression |>
    dplyr::left_join(metadata, by = c("cohort_id", "sample_id")) |>
    dplyr::select(
      .data$cohort_id, .data$sample_id, .data$gene_symbol,
      .data$expr_value, .data$expression_unit, .data$rank_percentile,
      .data$trial_id, .data$treatment_arm, .data$os_time_months,
      .data$os_event, .data$pfs_time_months, .data$pfs_event,
      .data$tumor_purity
    )

  if (nrow(analysis_long) != nrow(expression)) {
    stop("The expression-to-metadata join changed the number of rows.",
         call. = FALSE)
  }
  analysis_long
}
