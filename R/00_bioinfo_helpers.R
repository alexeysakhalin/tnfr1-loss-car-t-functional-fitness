# Shared validation and scoring helpers for the CheckMate analyses.

TCELL_SCORE_GENES <- c("CD2", "CD3D", "CD3E", "CD8A", "CD8B")

assert_columns <- function(x, required, object_name = deparse(substitute(x))) {
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop(object_name, " is missing required columns: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

safe_z <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  sx <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sx) || sx == 0) return(rep(0, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / sx)
}

validate_binary_event <- function(x, name) {
  values <- sort(unique(stats::na.omit(suppressWarnings(as.numeric(x)))))
  if (!identical(values, c(0, 1)) && !identical(values, 0) && !identical(values, 1)) {
    stop(name, " must contain only 0/1 (1 = observed event; 0 = censored); got: ",
         paste(values, collapse = ", "))
  }
  invisible(TRUE)
}

checkmate_nivolumab_metadata <- function(analysis_long) {
  required <- c(
    "cohort_id", "sample_id", "trial_id", "treatment_arm",
    "os_time_months", "os_event", "pfs_time_months", "pfs_event"
  )
  assert_columns(analysis_long, required, "analysis_long")

  meta <- analysis_long |>
    dplyr::filter(.data$cohort_id == "CHECKMATE_CCRCC") |>
    dplyr::transmute(
      sample_id = as.character(.data$sample_id),
      trial_id = as.character(.data$trial_id),
      treatment_arm = toupper(trimws(as.character(.data$treatment_arm))),
      os_time_months = suppressWarnings(as.numeric(.data$os_time_months)),
      os_event = suppressWarnings(as.numeric(.data$os_event)),
      pfs_time_months = suppressWarnings(as.numeric(.data$pfs_time_months)),
      pfs_event = suppressWarnings(as.numeric(.data$pfs_event))
    ) |>
    dplyr::distinct()

  duplicate_meta <- meta |>
    dplyr::count(.data$sample_id) |>
    dplyr::filter(.data$n != 1)
  if (nrow(duplicate_meta)) stop("Conflicting CheckMate metadata for one or more RNA samples.")
  if (nrow(meta) != 311L) stop("Expected 311 CheckMate RNA samples; found ", nrow(meta), ".")

  arm_counts <- table(meta$treatment_arm)
  expected_arms <- c(EVEROLIMUS = 130L, NIVOLUMAB = 181L)
  if (!identical(as.integer(arm_counts[names(expected_arms)]), unname(expected_arms))) {
    stop("Unexpected CheckMate arm counts. Expected EVEROLIMUS=130 and NIVOLUMAB=181; got ",
         paste(names(arm_counts), arm_counts, sep = "=", collapse = ", "), ".")
  }

  validate_binary_event(meta$os_event, "OS_CNSR/os_event")
  validate_binary_event(meta$pfs_event, "PFS_CNSR/pfs_event")
  if (any(!is.finite(meta$os_time_months)) || any(is.na(meta$os_event)) ||
      any(!is.finite(meta$pfs_time_months)) || any(is.na(meta$pfs_event))) {
    stop("CheckMate RNA samples require complete OS and PFS follow-up fields.")
  }
  if (sum(meta$os_event == 1, na.rm = TRUE) != 231L ||
      sum(meta$pfs_event == 1, na.rm = TRUE) != 276L) {
    stop("Unexpected all-arm event counts. Expected OS=231 and PFS=276 events.")
  }

  primary <- dplyr::filter(meta, .data$treatment_arm == "NIVOLUMAB")
  expected_trials <- c(`CM-009` = 16L, `CM-010` = 45L, `CM-025` = 120L)
  trial_counts <- table(primary$trial_id)
  if (!identical(as.integer(trial_counts[names(expected_trials)]), unname(expected_trials))) {
    stop("Unexpected nivolumab trial counts. Expected CM-009=16, CM-010=45, CM-025=120; got ",
         paste(names(trial_counts), trial_counts, sep = "=", collapse = ", "), ".")
  }
  if (sum(primary$os_event == 1, na.rm = TRUE) != 123L ||
      sum(primary$pfs_event == 1, na.rm = TRUE) != 159L) {
    stop("Unexpected nivolumab event counts. Expected OS=123 and PFS=159 events.")
  }
  primary
}

pick_marker_column <- function(nm, exact, pattern) {
  hit <- nm[tolower(nm) %in% tolower(exact)]
  if (length(hit)) return(hit[[1]])
  hit <- grep(pattern, nm, ignore.case = TRUE, value = TRUE)
  if (length(hit)) return(hit[[1]])
  NA_character_
}

read_cluster_gene_sets <- function(marker_file, top_n = 20L, excluded_genes = character()) {
  if (!file.exists(marker_file)) stop("Signature-definition file not found: ", marker_file)
  extension <- tolower(tools::file_ext(marker_file))
  markers <- if (extension == "csv") {
    readr::read_csv(marker_file, show_col_types = FALSE)
  } else if (extension %in% c("xlsx", "xls")) {
    readxl::read_excel(marker_file)
  } else {
    stop("Signature-definition file must be CSV or Excel: ", marker_file)
  }
  markers$input_order <- seq_len(nrow(markers))
  cluster_col <- pick_marker_column(names(markers),
                                    c("cluster", "cluster_id", "cluster_label"), "cluster")
  gene_col <- pick_marker_column(names(markers),
                                 c("gene", "gene_symbol", "symbol"), "gene|symbol")
  logfc_col <- pick_marker_column(names(markers),
                                  c("avg_log2FC", "avg_logFC", "log2FC", "logFC"), "log.*fc")
  padj_col <- pick_marker_column(names(markers),
                                 c("p_val_adj", "padj", "FDR"), "p.*adj|padj|fdr")
  if (is.na(cluster_col) || is.na(gene_col)) stop("Marker table lacks cluster/gene columns.")

  clean <- dplyr::tibble(
    cluster = paste0("C", gsub("^(Cluster_|cluster_|C)", "", as.character(markers[[cluster_col]]))),
    gene = as.character(markers[[gene_col]]),
    logfc = if (is.na(logfc_col)) NA_real_ else suppressWarnings(as.numeric(markers[[logfc_col]])),
    padj = if (is.na(padj_col)) NA_real_ else suppressWarnings(as.numeric(markers[[padj_col]])),
    input_order = markers$input_order
  ) |>
    dplyr::filter(
      .data$cluster %in% paste0("C", 0:9),
      !is.na(.data$gene), .data$gene != "",
      !.data$gene %in% excluded_genes
    ) |>
    dplyr::group_by(.data$cluster) |>
    dplyr::arrange(dplyr::desc(.data$logfc), .data$input_order, .by_group = TRUE) |>
    dplyr::distinct(.data$cluster, .data$gene, .keep_all = TRUE) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::ungroup()

  sets <- split(clean$gene, clean$cluster)
  missing <- setdiff(paste0("C", 0:9), names(sets))
  if (length(missing)) stop("Missing marker sets for: ", paste(missing, collapse = ", "))
  sets <- sets[paste0("C", 0:9)]
  if (any(lengths(sets) < 3L)) {
    stop("Each transferred state signature must contain at least three genes.")
  }
  sets
}

build_nivolumab_expression <- function(analysis_long, genes_needed) {
  primary_meta <- checkmate_nivolumab_metadata(analysis_long)
  analysis_long |>
    dplyr::filter(
      .data$cohort_id == "CHECKMATE_CCRCC",
      toupper(trimws(as.character(.data$treatment_arm))) == "NIVOLUMAB",
      .data$gene_symbol %in% genes_needed
    ) |>
    dplyr::transmute(
      sample_id = as.character(.data$sample_id),
      gene_symbol = as.character(.data$gene_symbol),
      expr_value = suppressWarnings(as.numeric(.data$expr_value))
    ) |>
    dplyr::group_by(.data$sample_id, .data$gene_symbol) |>
    dplyr::summarise(expr_value = mean(.data$expr_value, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = .data$gene_symbol, values_from = .data$expr_value) |>
    dplyr::left_join(primary_meta, by = "sample_id")
}

fit_trial_stratified_survival <- function(data, score_col, endpoint, signature) {
  time_col <- paste0(tolower(endpoint), "_time_months")
  event_col <- paste0(tolower(endpoint), "_event")
  dat <- data |>
    dplyr::transmute(
      score = suppressWarnings(as.numeric(.data[[score_col]])),
      time = suppressWarnings(as.numeric(.data[[time_col]])),
      event = suppressWarnings(as.numeric(.data[[event_col]])),
      trial_id = factor(.data$trial_id)
    ) |>
    dplyr::filter(is.finite(.data$score), is.finite(.data$time), !is.na(.data$event),
                  !is.na(.data$trial_id))
  validate_binary_event(dat$event, paste0(endpoint, " event"))
  cutoff <- stats::median(dat$score, na.rm = TRUE)
  dat$group <- factor(ifelse(dat$score > cutoff, "High", "Low"), levels = c("Low", "High"))
  dat$score_z <- safe_z(dat$score)

  split_fit <- survival::coxph(survival::Surv(time, event) ~ group + strata(trial_id), data = dat)
  continuous_fit <- survival::coxph(survival::Surv(time, event) ~ score_z + strata(trial_id), data = dat)
  # The KM/log-rank contrast is the unadjusted descriptive comparison shown in
  # the figure. Trial is handled in both Cox models through strata(trial_id).
  lr <- survival::survdiff(survival::Surv(time, event) ~ group, data = dat)
  split_tidy <- broom::tidy(split_fit, exponentiate = TRUE, conf.int = TRUE)
  continuous_tidy <- broom::tidy(continuous_fit, exponentiate = TRUE, conf.int = TRUE)
  # Match the independent Python validation, which uses the rank time
  # transformation for the Schoenfeld-residual proportional-hazards test.
  ph_table <- survival::cox.zph(split_fit, transform = "rank")$table
  ph_rows <- setdiff(rownames(ph_table), "GLOBAL")
  group_ph_row <- grep("^group(High)?$", ph_rows, value = TRUE)
  if (length(group_ph_row) != 1L) {
    stop(
      "Could not identify the group term in the proportional-hazards test: ",
      paste(ph_rows, collapse = ", ")
    )
  }
  ph_p <- unname(ph_table[group_ph_row, "p"])
  if (length(ph_p) != 1L || !is.finite(ph_p)) {
    stop("The proportional-hazards test returned a non-finite group p-value.")
  }

  result <- dplyr::tibble(
    endpoint = endpoint,
    signature = signature,
    n = nrow(dat),
    events = sum(dat$event == 1),
    n_low = sum(dat$group == "Low"),
    n_high = sum(dat$group == "High"),
    median_score = cutoff,
    logrank_p = stats::pchisq(lr$chisq, df = 1, lower.tail = FALSE),
    HR_high_vs_low = split_tidy$estimate[[1]],
    CI_low = split_tidy$conf.low[[1]],
    CI_high = split_tidy$conf.high[[1]],
    cox_p = split_tidy$p.value[[1]],
    continuous_HR_per_SD = continuous_tidy$estimate[[1]],
    continuous_CI_low = continuous_tidy$conf.low[[1]],
    continuous_CI_high = continuous_tidy$conf.high[[1]],
    continuous_cox_p = continuous_tidy$p.value[[1]],
    PH_test_p = ph_p,
    trial_stratified_cox = TRUE,
    analysis = "nivolumab_primary",
    arm = "NIVOLUMAB"
  )
  list(result = result, data = dat,
       km_fit = survival::survfit(survival::Surv(time, event) ~ group, data = dat))
}

save_ggsurvplot_pair <- function(plot_object, path, width = 7, height = 7.2,
                                 dpi = 600) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    png = path,
    tiff = paste0(tools::file_path_sans_ext(path), ".tiff")
  )
  for (format_name in names(paths)) {
    if (format_name == "png") {
      grDevices::png(
        filename = paths[[format_name]], width = width, height = height,
        units = "in", res = dpi, bg = "white"
      )
    } else {
      grDevices::tiff(
        filename = paths[[format_name]], width = width, height = height,
        units = "in", res = dpi, compression = "lzw", bg = "white"
      )
    }
    tryCatch(print(plot_object), finally = grDevices::dev.off())
  }
  invisible(paths)
}
