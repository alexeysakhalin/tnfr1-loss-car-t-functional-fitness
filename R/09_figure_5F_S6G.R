# Figure 5F and Supplementary Figure S6G: aggregate results from the
# exploratory CheckMate C6-context analysis.
#
# Patient-level expression is not required to render these panels. The local
# preparation workflow fits one model per eligible gene in 181
# nivolumab-treated tumors, excluding genes used to construct the C6 or bulk
# T-cell scores. Each model adjusts for the bulk T-cell score and trial. The
# complete aggregate model family is supplied here so BH correction is applied
# over every tested gene, not only the prespecified display genes.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
})

source(file.path("R", "00_bioinfo_helpers.R"))
source(file.path("R", "plot_style.R"))

global_model_file <- file.path(
  "data", "analysis", "checkmate_c6_global_gene_models.tsv.gz"
)
balance_file <- file.path(
  "data", "analysis", "checkmate_c6_group_balance.tsv"
)
curated_file <- file.path("resources", "Figure_5F_curated_gene_sets.csv")
output_dir <- file.path("results", "figure5F")
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

for (path in c(global_model_file, balance_file, curated_file)) {
  if (!file.exists(path)) stop("Required input is missing: ", path)
}

global_models <- read_tsv(
  global_model_file, show_col_types = FALSE, progress = FALSE
)
model_columns <- c(
  "gene", "n", "beta", "se", "ci_low", "ci_high", "p",
  "mean_high", "mean_low"
)
assert_columns(global_models, model_columns, "C6 global gene-model table")
global_models <- global_models |>
  transmute(
    gene = trimws(as.character(.data$gene)),
    n = suppressWarnings(as.integer(.data$n)),
    beta = suppressWarnings(as.numeric(.data$beta)),
    se = suppressWarnings(as.numeric(.data$se)),
    ci_low = suppressWarnings(as.numeric(.data$ci_low)),
    ci_high = suppressWarnings(as.numeric(.data$ci_high)),
    p = suppressWarnings(as.numeric(.data$p)),
    mean_high = suppressWarnings(as.numeric(.data$mean_high)),
    mean_low = suppressWarnings(as.numeric(.data$mean_low))
  )

if (nrow(global_models) < 1000L) {
  stop("The aggregate input is too small to represent the tested gene family.")
}
if (any(is.na(global_models$gene) | global_models$gene == "") ||
    anyDuplicated(global_models$gene)) {
  stop("The aggregate input must contain one uniquely named row per tested gene.")
}
numeric_columns <- setdiff(model_columns, "gene")
if (any(!vapply(global_models[numeric_columns], function(x) all(is.finite(x)), logical(1)))) {
  stop("The aggregate gene-model table contains non-finite model statistics.")
}
if (any(global_models$p < 0 | global_models$p > 1) ||
    any(global_models$n < 163L | global_models$n > 181L) ||
    any(global_models$se < 0) ||
    any(global_models$ci_low > global_models$ci_high) ||
    any(global_models$beta < global_models$ci_low |
        global_models$beta > global_models$ci_high)) {
  stop("The aggregate gene-model table failed range checks.")
}

# The adjustment is deliberately calculated only after the complete aggregate
# family has passed validation.
global_models <- global_models |>
  mutate(
    BH_p = p.adjust(.data$p, method = "BH"),
    neglog10_p = -log10(pmax(.data$p, .Machine$double.xmin)),
    neglog10_BH_p = -log10(pmax(.data$BH_p, .Machine$double.xmin)),
    direction = if_else(
      .data$beta >= 0,
      "higher in adjusted C6-high", "higher in adjusted C6-low"
    )
  ) |>
  arrange(.data$p)
write_tsv(
  global_models,
  file.path(table_dir, "Figure_5F_global_gene_models_with_BH.tsv.gz")
)

curated <- read_csv(curated_file, show_col_types = FALSE) |>
  transmute(
    gene = trimws(as.character(.data$gene)),
    block = trimws(as.character(.data$block))
  ) |>
  distinct()
assert_columns(curated, c("gene", "block"), "Figure 5F display set")
if (anyDuplicated(curated$gene)) {
  stop("Each gene may appear only once in the Figure 5F display set.")
}
missing_curated <- setdiff(curated$gene, global_models$gene)
if (length(missing_curated)) {
  stop(
    "Display genes absent from the complete model family: ",
    paste(missing_curated, collapse = ", ")
  )
}
curated_statistics <- curated |>
  left_join(global_models, by = "gene") |>
  mutate(selection = "Prespecified display set; statistics from global model")
write_tsv(
  curated_statistics,
  file.path(table_dir, "Figure_5F_curated_gene_statistics.tsv")
)

block_colors <- c(
  "Apoptosis / stress" = "#D62728",
  "Survival / adaptation" = "#3182BD"
)
unknown_blocks <- setdiff(unique(curated_statistics$block), names(block_colors))
if (length(unknown_blocks)) {
  stop("Unrecognized Figure 5F display blocks: ", paste(unknown_blocks, collapse = ", "))
}

volcano_plot <- ggplot(global_models, aes(.data$beta, .data$neglog10_BH_p)) +
  geom_point(color = "grey78", alpha = 0.35, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey45") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
  geom_point(
    data = curated_statistics,
    aes(color = .data$block), size = 3.1
  ) +
  ggrepel::geom_text_repel(
    data = curated_statistics,
    aes(label = .data$gene, color = .data$block),
    size = 3.2, max.overlaps = Inf, min.segment.length = 0,
    box.padding = 0.35, show.legend = FALSE
  ) +
  scale_color_manual(values = block_colors) +
  labs(
    title = "Exploratory C6-associated bulk expression context",
    subtitle = paste0(
      "Nivolumab-treated ccRCC (n=181); adjusted for bulk T-cell score and trial; ",
      "BH correction across all fitted genes"
    ),
    x = "Adjusted coefficient: C6-high versus C6-low",
    y = expression(-log[10]("BH-adjusted p")), color = NULL
  ) +
  theme_pub() +
  theme(legend.position = "top")

save_ggplot_pair(
  volcano_plot,
  file.path(figure_dir, "Figure_5F_C6_adjusted_gene_level_volcano.png"),
  width = 9, height = 6.5
)

# S6G uses the two-group aggregate summary from the same 181-sample model.
balance <- read_tsv(balance_file, show_col_types = FALSE, progress = FALSE)
balance_columns <- c(
  "c6_group", "n", "tcell_score_mean", "tcell_score_sd", "welch_p"
)
assert_columns(balance, balance_columns, "C6 group-balance table")
balance <- balance |>
  transmute(
    c6_group = trimws(as.character(.data$c6_group)),
    n = suppressWarnings(as.integer(.data$n)),
    tcell_score_mean = suppressWarnings(as.numeric(.data$tcell_score_mean)),
    tcell_score_sd = suppressWarnings(as.numeric(.data$tcell_score_sd)),
    welch_p = suppressWarnings(as.numeric(.data$welch_p))
  )

if (nrow(balance) != 2L || anyDuplicated(balance$c6_group) ||
    !setequal(balance$c6_group, c("Low", "High"))) {
  stop("The balance table must contain exactly one Low and one High row.")
}
if (sum(balance$n) != 181L || any(balance$n <= 0L)) {
  stop("The Low and High balance groups must contain 181 samples in total.")
}
if (any(!is.finite(balance$tcell_score_mean)) ||
    any(!is.finite(balance$tcell_score_sd)) ||
    any(balance$tcell_score_sd < 0)) {
  stop("The balance table contains invalid score summaries.")
}
welch_values <- unique(balance$welch_p)
if (length(welch_values) != 1L || !is.finite(welch_values) ||
    welch_values < 0 || welch_values > 1) {
  stop("The balance table must contain one valid Welch-test p value.")
}
balance <- balance |>
  mutate(
    c6_group = factor(.data$c6_group, levels = c("Low", "High")),
    lower = .data$tcell_score_mean - .data$tcell_score_sd,
    upper = .data$tcell_score_mean + .data$tcell_score_sd
  )
write_tsv(balance, file.path(table_dir, "Supplementary_Figure_S6G_balance.tsv"))

balance_plot <- ggplot(
  balance,
  aes(.data$c6_group, .data$tcell_score_mean, color = .data$c6_group)
) +
  geom_hline(yintercept = 0, color = "grey75", linewidth = 0.5) +
  geom_errorbar(
    aes(ymin = .data$lower, ymax = .data$upper),
    width = 0.12, linewidth = 0.8
  ) +
  geom_point(size = 4) +
  geom_text(
    aes(label = paste0("n = ", .data$n), y = .data$upper),
    nudge_y = 0.12, show.legend = FALSE
  ) +
  scale_color_manual(values = c(Low = "grey55", High = "#29AFC4")) +
  labs(
    title = "Residualized C6 group balance",
    subtitle = sprintf(
      "Mean +/- SD; Welch p=%.3g; descriptive balance check", welch_values
    ),
    x = "Adjusted C6 group", y = "Bulk T-cell expression score (z)"
  ) +
  theme_pub() +
  theme(legend.position = "none")

save_ggplot_pair(
  balance_plot,
  file.path(figure_dir, "Supplementary_Figure_S6G_C6_group_balance.png"),
  width = 5.5, height = 5
)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
