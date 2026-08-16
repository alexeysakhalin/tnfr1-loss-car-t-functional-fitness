# Figures 5E/S6D/S6E: transferred in-vitro T-cell-state expression signatures.
# Primary population: NIVOLUMAB only (n=181). Cox models are stratified by
# trial. BH correction is applied over the complete family of ten signatures
# (C0-C9), separately for endpoint and test. CNSR coding is 1=event, 0=censored.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(survival)
  library(survminer)
  library(broom)
})

source(file.path("R", "00_bioinfo_helpers.R"))
source(file.path("R", "00_load_analysis_tables.R"))
source(file.path("R", "plot_style.R"))

marker_file <- file.path("resources", "CAR_T_state_signatures.csv")
out_dir <- "figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

analysis_long <- load_analysis_tables()
cluster_sets <- read_cluster_gene_sets(marker_file, top_n = 20L)
genes_needed <- unique(unlist(cluster_sets, use.names = FALSE))
dat <- build_nivolumab_expression(analysis_long, genes_needed)
if (nrow(dat) != 181L || anyDuplicated(dat$sample_id)) {
  stop("Figures 5E/S6D/S6E require 181 unique nivolumab RNA samples.")
}

coverage <- purrr::imap_dfr(cluster_sets, function(genes, cluster) {
  present <- intersect(genes, names(dat))
  tibble(
    analysis_arm = "NIVOLUMAB", signature = cluster,
    n_genes_defined = length(genes), n_genes_present = length(present),
    genes_present = paste(present, collapse = ";")
  )
})
if (any(coverage$n_genes_present < 3L)) {
  stop("At least one signature has fewer than three genes in bulk expression.")
}
write.csv(coverage, file.path(out_dir, "CheckMate_signature_coverage_nivolumab.csv"),
          row.names = FALSE)

all_present <- intersect(genes_needed, names(dat))
if (any(!is.finite(as.matrix(dat[all_present])))) {
  stop("Transferred-state scores require complete expression for present genes.")
}
dat[all_present] <- lapply(dat[all_present], safe_z)
for (cluster in names(cluster_sets)) {
  present <- intersect(cluster_sets[[cluster]], names(dat))
  dat[[paste0(cluster, "_score")]] <- rowMeans(dat[present], na.rm = FALSE)
}

fits <- list()
results <- list()
for (endpoint in c("OS", "PFS")) {
  for (cluster in paste0("C", 0:9)) {
    key <- paste(endpoint, cluster, sep = "_")
    fits[[key]] <- fit_trial_stratified_survival(
      dat, paste0(cluster, "_score"), endpoint, cluster
    )
    results[[key]] <- fits[[key]]$result
  }
}
results <- bind_rows(results) |>
  group_by(endpoint) |>
  mutate(
    logrank_BH_10_states = p.adjust(logrank_p, method = "BH"),
    cox_BH_10_states = p.adjust(cox_p, method = "BH"),
    continuous_cox_BH_10_states = p.adjust(continuous_cox_p, method = "BH")
  ) |>
  ungroup() |>
  arrange(endpoint, signature)

stopifnot(all(table(results$endpoint) == 10L))
write.csv(results,
          file.path(out_dir, "CheckMate_survival_results_nivolumab_primary.csv"),
          row.names = FALSE)

cluster_cols <- c(
  C0 = "#F08A80", C1 = "#D98C00", C2 = "#B8A500", C3 = "#6CB400",
  C4 = "#20B95A", C5 = "#16B8A8", C6 = "#29AFC4", C7 = "#2A96E6",
  C8 = "#A27AE8", C9 = "#D865D8"
)

format_p <- function(x) {
  if (is.na(x)) return("NA")
  if (x < 0.001) return("< 0.001")
  sprintf("%.3f", x)
}

save_km <- function(endpoint, cluster, filename) {
  key <- paste(endpoint, cluster, sep = "_")
  fit <- fits[[key]]
  r <- results |>
    filter(
      .data$endpoint == .env$endpoint,
      .data$signature == .env$cluster
    )
  if (nrow(r) != 1L) {
    stop(
      sprintf(
        "Expected exactly one survival result for endpoint=%s and signature=%s; found %d.",
        endpoint, cluster, nrow(r)
      ),
      call. = FALSE
    )
  }
  endpoint_label <- if (endpoint == "OS") {
    "Overall survival"
  } else {
    "Progression-free survival"
  }
  subtitle <- sprintf(
    paste0(
      "Nivolumab only (n = %d)\n",
      "Trial-stratified Cox (median split): HR %.2f (95%% CI %.2f-%.2f)\n",
      "Cox p = %s; BH p = %s\n",
      "Median-split log-rank (descriptive): p = %s; BH p = %s"
    ),
    nrow(fit$data), r$HR_high_vs_low, r$CI_low, r$CI_high,
    format_p(r$cox_p), format_p(r$cox_BH_10_states),
    format_p(r$logrank_p), format_p(r$logrank_BH_10_states)
  )
  p <- survminer::ggsurvplot(
    fit$km_fit, data = fit$data, risk.table = TRUE, censor = TRUE,
    palette = c("grey65", unname(cluster_cols[[cluster]])),
    legend.title = paste(cluster, "score"), legend.labs = c("Low", "High"),
    title = paste(endpoint_label, "by", cluster, "transferred T-cell-state signature"),
    subtitle = subtitle, xlab = "Time (months)", ylab = "Survival probability"
  )
  p$plot <- p$plot +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 17, face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 11.5, lineheight = 1.12, hjust = 0),
      plot.margin = ggplot2::margin(t = 12, r = 20, b = 8, l = 12)
    )
  save_ggsurvplot_pair(p, file.path(out_dir, filename))
}

save_km("OS", "C6", "Figure_5E_KM_ccRCC_C6.png")
save_km("OS", "C1", "Supplementary_Figure_S6D_C1.png")
save_km("OS", "C0", "Supplementary_Figure_S6D_C0.png")
save_km("PFS", "C6", "Supplementary_Figure_S6E_PFS_C6.png")
save_km("PFS", "C0", "Supplementary_Figure_S6E_PFS_C0.png")
writeLines(
  capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_R08.txt")
)
