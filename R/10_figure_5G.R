# Figure 5G: B2M/ICAM1 association with transferred in-vitro T-cell-state scores.
# Scores are bulk expression signatures, not cell fractions. Primary analysis
# uses nivolumab-treated CheckMate samples only (n=181) and adjusts for total
# T-cell score plus trial. Because individual ICAM1, individual B2M and their
# combined score were evaluated, BH correction covers all 30 fitted models.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(broom)
})

source(file.path("R", "00_bioinfo_helpers.R"))
source(file.path("R", "00_load_analysis_tables.R"))

marker_file <- file.path("resources", "CAR_T_state_signatures.csv")
out_dir <- "figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tcell_genes <- TCELL_SCORE_GENES
predictor_genes <- c("ICAM1", "B2M")
clusters <- paste0("C", 0:9)
excluded <- unique(c(tcell_genes, predictor_genes))
cluster_sets <- read_cluster_gene_sets(marker_file, top_n = 20L,
                                       excluded_genes = excluded)

analysis_long <- load_analysis_tables()
genes_needed <- unique(c(unlist(cluster_sets, use.names = FALSE),
                         tcell_genes, predictor_genes))
dat <- build_nivolumab_expression(analysis_long, genes_needed)
if (nrow(dat) != 181L) stop("Figure 5G primary analysis must contain 181 samples.")

tcell_present <- intersect(tcell_genes, names(dat))
if (!setequal(tcell_present, tcell_genes) ||
    !all(predictor_genes %in% names(dat))) {
  stop("Insufficient T-cell genes or missing ICAM1/B2M.")
}
signature_present <- intersect(unique(unlist(cluster_sets)), names(dat))
all_model_genes <- unique(c(tcell_present, predictor_genes, signature_present))
if (any(!is.finite(as.matrix(dat[all_model_genes])))) {
  stop("Figure 5G requires complete expression for all genes used in a model.")
}
dat[all_model_genes] <- lapply(dat[all_model_genes], safe_z)

coverage <- imap_dfr(cluster_sets, function(genes, cluster) {
  present <- intersect(genes, names(dat))
  tibble(
    cluster = cluster, n_signature_genes_defined = length(genes),
    n_signature_genes_present = length(present),
    genes_present = paste(present, collapse = ";")
  )
})
if (any(coverage$n_signature_genes_present < 3L)) {
  stop("At least one Figure 5G signature has fewer than three genes.")
}
write.csv(coverage, file.path(out_dir, "Figure_5G_signature_definitions.csv"),
          row.names = FALSE)

for (cluster in clusters) {
  present <- intersect(cluster_sets[[cluster]], names(dat))
  dat[[paste0(cluster, "_state")]] <- safe_z(rowMeans(dat[present], na.rm = FALSE))
}
dat <- dat |>
  mutate(
    TcellScore_z = safe_z(rowMeans(across(all_of(tcell_present)), na.rm = FALSE)),
    ICAM1_z = safe_z(.data$ICAM1),
    B2M_z = safe_z(.data$B2M),
    combined_B2M_ICAM1_z = safe_z((.data$ICAM1_z + .data$B2M_z) / 2),
    trial_id = factor(.data$trial_id)
  )

fit_signature <- function(cluster, predictor) {
  outcome <- paste0(cluster, "_state")
  tmp <- dat |>
    select(all_of(c(outcome, predictor, "TcellScore_z", "trial_id"))) |>
    filter(if_all(-trial_id, is.finite), !is.na(.data$trial_id))
  fit <- lm(as.formula(paste(outcome, "~", predictor,
                             "+ TcellScore_z + factor(trial_id)")), data = tmp)
  broom::tidy(fit, conf.int = TRUE) |>
    filter(.data$term == predictor) |>
    transmute(
      n = nrow(tmp), beta = .data$estimate, se = .data$std.error,
      ci_low = .data$conf.low, ci_high = .data$conf.high, p = .data$p.value,
      df_resid = stats::df.residual(fit), cluster = cluster,
      predictor = predictor, covariates = "TcellScore_z + trial",
      n_signature_genes = coverage$n_signature_genes_present[coverage$cluster == cluster]
    )
}

predictors <- c("ICAM1_z", "B2M_z", "combined_B2M_ICAM1_z")
results <- map_dfr(predictors, function(predictor) {
  map_dfr(clusters, fit_signature, predictor = predictor)
}) |>
  mutate(BH_p_30_models = p.adjust(.data$p, method = "BH")) |>
  arrange(.data$predictor, .data$cluster)

if (any(results$n != 181L) || nrow(results) != 30L ||
    any(table(results$predictor) != 10L)) {
  stop("Figure 5G model family is incomplete or sample counts differ from 181.")
}
write.csv(results,
          file.path(out_dir, "Figure_5G_signature_association_results.csv"),
          row.names = FALSE)

score_cols <- paste0(clusters, "_state")
write.csv(cor(dat[score_cols], use = "pairwise.complete.obs"),
          file.path(out_dir, "Figure_5G_signature_score_correlations.csv"))

cluster_cols <- c(
  C0 = "#F08A80", C1 = "#D98C00", C2 = "#B8A500", C3 = "#6CB400",
  C4 = "#20B95A", C5 = "#16B8A8", C6 = "#29AFC4", C7 = "#2A96E6",
  C8 = "#A27AE8", C9 = "#D865D8"
)
plot_df <- results |>
  filter(.data$predictor == "combined_B2M_ICAM1_z") |>
  mutate(cluster = factor(.data$cluster, levels = clusters))

p <- ggplot(plot_df, aes(.data$cluster, .data$beta, fill = .data$cluster)) +
  geom_col(width = 0.78) +
  geom_errorbar(aes(ymin = .data$ci_low, ymax = .data$ci_high), width = 0.18) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(y = .data$ci_high + 0.045,
                label = ifelse(.data$BH_p_30_models < 0.001, "BH p<0.001",
                               sprintf("BH p=%.3f", .data$BH_p_30_models))), size = 3) +
  scale_fill_manual(values = cluster_cols, drop = FALSE) +
  labs(
    title = "B2M+ICAM1 association with transferred T-cell-state scores",
    subtitle = paste0("Nivolumab-treated ccRCC (n=181); adjusted for T-cell score and trial; ",
                      "BH across 30 models"),
    x = "Tumor-co-culture-derived state score", y = "Adjusted standardized beta"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5), plot.margin = margin(10, 20, 10, 10))

ggsave(file.path(out_dir, "Figure_5G_B2M_ICAM1_signature_associations.png"),
       p, width = 8.5, height = 5.5, dpi = 600, bg = "white")
