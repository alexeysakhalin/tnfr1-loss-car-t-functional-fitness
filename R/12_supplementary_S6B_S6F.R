# Supplementary Figures S6B and S6F: exploratory bulk-expression correlations.
#
# All three score definitions use the same gene members in all four cohorts.
# CheckMate is restricted to the 181 nivolumab-treated tumors. The 12 Spearman
# tests (three score pairs by four cohorts) form one BH-corrected family.
# Scores are within-cohort standardized bulk-expression summaries, not cell
# fractions.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(ggplot2)
})

source(file.path("R", "00_bioinfo_helpers.R"))
source(file.path("R", "00_load_analysis_tables.R"))
source(file.path("R", "plot_style.R"))

out_dir <- file.path("results", "supplementary_S6")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

analysis_long <- load_analysis_tables()
checkmate_nivolumab_metadata(analysis_long)

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

defined_sets <- list(
  Tcell = TCELL_SCORE_GENES,
  TNF_IFNG = c("IFNG", "IRF1", "STAT1", "CXCL9", "CXCL10", "ICAM1", "TNFRSF1A"),
  NFkB = c("NFKBIA", "TNFAIP3", "RELB", "BIRC3", "CXCL2", "CXCL3",
           "JUNB", "FOS", "EGR1", "TRAF1", "NFKB2", "IER3",
           "PTGS2", "IL6", "CCL2", "CCL20", "TNIP1")
)
all_defined <- unique(c(unlist(defined_sets), "ICAM1"))

selected_expression <- analysis_long |>
  filter(
    .data$cohort_id %in% cohort_order,
    .data$gene_symbol %in% all_defined,
    .data$cohort_id != "CHECKMATE_CCRCC" |
      toupper(trimws(as.character(.data$treatment_arm))) == "NIVOLUMAB"
  ) |>
  transmute(
    cohort_id = as.character(.data$cohort_id),
    sample_id = as.character(.data$sample_id),
    gene = as.character(.data$gene_symbol),
    expression = suppressWarnings(as.numeric(.data$expr_value))
  ) |>
  filter(is.finite(.data$expression)) |>
  group_by(.data$cohort_id, .data$sample_id, .data$gene) |>
  summarise(expression = mean(.data$expression), .groups = "drop")

n_checkmate <- selected_expression |>
  filter(.data$cohort_id == "CHECKMATE_CCRCC") |>
  summarise(n = n_distinct(.data$sample_id)) |>
  pull(n)
if (length(n_checkmate) != 1L || n_checkmate != 181L) {
  stop("S6B/S6F CheckMate input must contain 181 nivolumab tumors.")
}

gene_by_cohort <- selected_expression |>
  distinct(.data$cohort_id, .data$gene) |>
  count(.data$gene, name = "n_cohorts")
common_genes <- gene_by_cohort |>
  filter(.data$n_cohorts == length(cohort_order)) |>
  pull(.data$gene)
common_sets <- map(defined_sets, intersect, y = common_genes)
minimums <- c(Tcell = 4L, TNF_IFNG = 4L, NFkB = 8L)
if (any(lengths(common_sets) < minimums[names(common_sets)]) || !"ICAM1" %in% common_genes) {
  stop("Insufficient common score-gene coverage across the four cohorts.")
}

manifest <- imap_dfr(common_sets, function(genes, score) {
  tibble(
    score = score,
    n_defined = length(defined_sets[[score]]),
    n_used_in_every_cohort = length(genes),
    genes = paste(genes, collapse = ";")
  )
})
write_csv(manifest, file.path(out_dir, "Supplementary_Figure_S6_score_manifest.csv"))

make_cohort_scores <- function(cohort) {
  wide <- selected_expression |>
    filter(.data$cohort_id == cohort, .data$gene %in% common_genes) |>
    select(.data$sample_id, .data$gene, .data$expression) |>
    pivot_wider(names_from = .data$gene, values_from = .data$expression)
  genes_to_scale <- intersect(common_genes, names(wide))
  wide <- wide |>
    mutate(across(all_of(genes_to_scale), safe_z)) |>
    mutate(
      Tcell_score = rowMeans(across(all_of(common_sets$Tcell)), na.rm = FALSE),
      TNF_IFNG_score = rowMeans(across(all_of(common_sets$TNF_IFNG)), na.rm = FALSE),
      NFkB_score = rowMeans(across(all_of(common_sets$NFkB)), na.rm = FALSE),
      ICAM1_z = .data$ICAM1,
      cohort_id = cohort,
      cohort_label = unname(cohort_labels[[cohort]])
    ) |>
    filter(if_all(c("Tcell_score", "TNF_IFNG_score", "NFkB_score", "ICAM1_z"),
                  is.finite))
  wide
}

score_df <- map_dfr(cohort_order, make_cohort_scores) |>
  mutate(cohort_label = factor(.data$cohort_label,
                               levels = unname(cohort_labels[cohort_order])))
checkmate_score_n <- score_df |>
  filter(.data$cohort_id == "CHECKMATE_CCRCC") |>
  summarise(n = n_distinct(.data$sample_id)) |>
  pull(.data$n)
if (!identical(checkmate_score_n, 181L)) {
  stop("S6B/S6F scores must be complete for all 181 nivolumab tumors.")
}

pairs <- tribble(
  ~comparison, ~x, ~y, ~x_label, ~y_label,
  "T-cell vs TNF/IFNG", "Tcell_score", "TNF_IFNG_score",
  "Bulk T-cell score", "TNF/IFNG score",
  "ICAM1 vs NF-kB", "ICAM1_z", "NFkB_score",
  "ICAM1 expression (z)", "NF-kB score",
  "TNF/IFNG vs NF-kB", "TNF_IFNG_score", "NFkB_score",
  "TNF/IFNG score", "NF-kB score"
)

cor_results <- crossing(
  cohort_id = cohort_order,
  comparison = pairs$comparison
) |>
  left_join(pairs, by = "comparison") |>
  mutate(result = map2(.data$cohort_id, .data$comparison, function(cohort, comparison) {
    spec <- pairs |> filter(.data$comparison == .env$comparison)
    dat <- score_df |> filter(.data$cohort_id == cohort)
    test <- suppressWarnings(cor.test(
      dat[[spec$x[[1]]]], dat[[spec$y[[1]]]],
      method = "spearman", exact = FALSE
    ))
    tibble(
      n = nrow(dat), rho = unname(test$estimate), p = test$p.value,
      cohort_label = unname(cohort_labels[[cohort]])
    )
  })) |>
  unnest(.data$result) |>
  mutate(BH_p_12_tests = p.adjust(.data$p, method = "BH")) |>
  arrange(match(.data$comparison, pairs$comparison), match(.data$cohort_id, cohort_order))

if (nrow(cor_results) != 12L) stop("Expected 12 correlation tests.")
write_csv(cor_results, file.path(out_dir, "Supplementary_Figure_S6_correlations.csv"))

cor_results <- cor_results |>
  mutate(
    cohort_label = factor(.data$cohort_label,
                          levels = unname(cohort_labels[cohort_order])),
    comparison = factor(.data$comparison, levels = pairs$comparison)
  )

p_bars <- ggplot(cor_results, aes(.data$cohort_label, .data$rho,
                                  fill = .data$cohort_label)) +
  geom_col(width = 0.72) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = sprintf("n=%d\nBH p=%.3g", .data$n,
                                .data$BH_p_12_tests),
                y = if_else(.data$rho >= 0, .data$rho + 0.08, .data$rho - 0.08)),
            size = 3) +
  facet_wrap(~ comparison, nrow = 1) +
  scale_fill_manual(values = cohort_colors) +
  coord_cartesian(ylim = c(-1.25, 1.25), clip = "off") +
  labs(
    title = "Exploratory within-cohort bulk-expression correlations",
    subtitle = "BH adjustment across 12 prespecified Spearman tests",
    x = NULL, y = "Spearman rho"
  ) +
  theme_pub() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(
  file.path(out_dir, "Supplementary_Figure_S6B_bulk_score_correlations.png"),
  p_bars, width = 13, height = 5.2, dpi = 600, bg = "white"
)

ccrcc_long <- pairs |>
  mutate(data = pmap(list(.data$comparison, .data$x, .data$y,
                          .data$x_label, .data$y_label),
                     function(comparison, x, y, x_label, y_label) {
    score_df |>
      filter(.data$cohort_id == "CHECKMATE_CCRCC") |>
      transmute(
        comparison = .env$comparison, x_value = .data[[x]], y_value = .data[[y]],
        x_label = .env$x_label, y_label = .env$y_label
      )
  })) |>
  select(.data$data) |>
  unnest(.data$data) |>
  mutate(comparison = factor(.data$comparison, levels = pairs$comparison))

ccrcc_labels <- cor_results |>
  filter(.data$cohort_id == "CHECKMATE_CCRCC") |>
  transmute(
    comparison = .data$comparison,
    label = sprintf("rho=%.2f; nominal p=%.3g; BH p=%.3g; n=%d",
                    .data$rho, .data$p, .data$BH_p_12_tests, .data$n)
  )

p_scatter <- ggplot(ccrcc_long, aes(.data$x_value, .data$y_value)) +
  geom_point(alpha = 0.55, size = 1.7, color = "#D62728") +
  geom_smooth(
    method = "loess", formula = y ~ x, span = 0.75,
    se = TRUE, color = "black", linewidth = 0.7
  ) +
  facet_wrap(~ comparison, scales = "free", nrow = 1) +
  geom_label(
    data = ccrcc_labels,
    aes(x = -Inf, y = Inf, label = .data$label),
    inherit.aes = FALSE, hjust = -0.02, vjust = 1.1, size = 3
  ) +
  labs(
    title = "Nivolumab-treated CheckMate ccRCC bulk-score correlations",
    subtitle = paste0(
      "Exploratory Spearman analyses; LOESS curves are descriptive; ",
      "panels use n=181 complete-score tumors"
    ),
    x = "First score (within-cohort standardized)",
    y = "Second score (within-cohort standardized)"
  ) +
  theme_pub()

ggsave(
  file.path(out_dir, "Supplementary_Figure_S6F_CheckMate_correlations.png"),
  p_scatter, width = 13, height = 5.2, dpi = 600, bg = "white"
)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
