# Supplementary Figure S1B: release-locked DepMap RIPK3/NLRP3 analysis.
#
# This script intentionally stops unless the exact DepMap release, DOI,
# download date and SHA-256 checksums are provided. It retains one official
# default expression entry per human tumor-derived cell-line model. The same
# fixed cutoff (0.5 log2(TPM+1)) is used for statistics and plotted dividers.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(digest)
  library(jsonlite)
})

data_dir <- file.path("data", "depmap")
out_dir <- file.path("results", "supplementary_S1B")
expr_file <- file.path(data_dir, "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
model_file <- file.path(data_dir, "Model.csv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- c(
  release = Sys.getenv("DEPMAP_RELEASE"),
  release_doi = Sys.getenv("DEPMAP_RELEASE_DOI"),
  download_date = Sys.getenv("DEPMAP_DOWNLOAD_DATE"),
  expression_sha256 = Sys.getenv("DEPMAP_EXPRESSION_SHA256"),
  model_sha256 = Sys.getenv("DEPMAP_MODEL_SHA256")
)
if (any(!nzchar(metadata))) {
  stop(
    "Set DEPMAP_RELEASE, DEPMAP_RELEASE_DOI, DEPMAP_DOWNLOAD_DATE, ",
    "DEPMAP_EXPRESSION_SHA256 and DEPMAP_MODEL_SHA256 before running."
  )
}
if (!grepl("^10\\.[0-9]{4,9}/", metadata[["release_doi"]])) {
  stop("DEPMAP_RELEASE_DOI is not DOI-shaped.")
}
if (!grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", metadata[["download_date"]])) {
  stop("DEPMAP_DOWNLOAD_DATE must use YYYY-MM-DD.")
}

assert_checksum <- function(path, expected) {
  if (!file.exists(path)) stop("Missing input: ", path)
  if (!grepl("^[0-9A-Fa-f]{64}$", expected)) stop("Invalid SHA-256 for ", path)
  observed <- digest::digest(path, algo = "sha256", file = TRUE)
  if (tolower(observed) != tolower(expected)) {
    stop("SHA-256 mismatch for ", path, "; observed ", observed)
  }
  observed
}
observed_expression_sha256 <- assert_checksum(expr_file, metadata[["expression_sha256"]])
observed_model_sha256 <- assert_checksum(model_file, metadata[["model_sha256"]])

expr_header <- names(fread(expr_file, nrows = 0, showProgress = FALSE))
model_header <- names(fread(model_file, nrows = 0, showProgress = FALSE))
required_model <- c("ModelID", "ModelType", "TissueOrigin", "OncotreePrimaryDisease")
if (length(setdiff(required_model, model_header))) {
  stop("Model.csv lacks required human tumor cell-line fields: ",
       paste(setdiff(required_model, model_header), collapse = ", "))
}
default_candidates <- c(
  "IsDefaultEntryForModel", "IsDefaultEntry", "isDefaultEntryForModel"
)
default_flag <- intersect(default_candidates, expr_header)[1]
if (is.na(default_flag) || !"ModelID" %in% expr_header) {
  stop("Expression file lacks ModelID or IsDefaultEntryForModel.")
}

find_gene_column <- function(symbol) {
  hits <- grep(paste0("^", symbol, "(\\s|$|\\()"), expr_header, value = TRUE)
  if (length(hits) != 1L) {
    stop("Expected one expression column for ", symbol, "; found ", length(hits), ".")
  }
  hits[[1]]
}
gene_columns <- setNames(
  vapply(c("RIPK3", "NLRP3"), find_gene_column, character(1)),
  c("RIPK3", "NLRP3")
)

eligible_models <- fread(
  model_file, select = required_model, showProgress = FALSE
) |>
  as_tibble() |>
  transmute(
    ModelID = as.character(.data$ModelID),
    ModelType = trimws(as.character(.data$ModelType)),
    TissueOrigin = trimws(as.character(.data$TissueOrigin)),
    OncotreePrimaryDisease = trimws(as.character(.data$OncotreePrimaryDisease))
  ) |>
  filter(
    toupper(.data$ModelType) == "CELL LINE",
    toupper(.data$TissueOrigin) == "HUMAN",
    !is.na(.data$OncotreePrimaryDisease), .data$OncotreePrimaryDisease != ""
  ) |>
  distinct(.data$ModelID, .keep_all = TRUE)
if (!nrow(eligible_models)) stop("No eligible human tumor cell-line models found.")

expression <- fread(
  expr_file,
  select = c("ModelID", default_flag, unname(gene_columns)),
  showProgress = FALSE
) |>
  as_tibble() |>
  transmute(
    ModelID = as.character(.data$ModelID),
    default_entry = toupper(trimws(as.character(.data[[default_flag]]))),
    RIPK3 = suppressWarnings(as.numeric(.data[[gene_columns[["RIPK3"]]]])),
    NLRP3 = suppressWarnings(as.numeric(.data[[gene_columns[["NLRP3"]]]]))
  ) |>
  filter(.data$default_entry %in% c("YES", "TRUE", "1")) |>
  select(-default_entry)
if (anyDuplicated(expression$ModelID)) {
  stop("Multiple default expression entries remain for at least one ModelID.")
}

analysis <- eligible_models |>
  inner_join(expression, by = "ModelID") |>
  filter(is.finite(.data$RIPK3), is.finite(.data$NLRP3))
if (nrow(analysis) < 100L) stop("Unexpectedly few eligible DepMap models.")

cutoff <- 0.5
analysis <- analysis |>
  mutate(
    RIPK3_low = .data$RIPK3 < cutoff,
    NLRP3_low = .data$NLRP3 < cutoff,
    quadrant = case_when(
      .data$RIPK3_low & .data$NLRP3_low ~ "Both low",
      .data$RIPK3_low & !.data$NLRP3_low ~ "RIPK3 low only",
      !.data$RIPK3_low & .data$NLRP3_low ~ "NLRP3 low only",
      TRUE ~ "Neither low"
    )
  )

summary_table <- bind_rows(
  tibble(
    metric = c("RIPK3 low", "NLRP3 low"),
    n = c(sum(analysis$RIPK3_low), sum(analysis$NLRP3_low)),
    denominator = nrow(analysis)
  ),
  analysis |>
    count(.data$quadrant, name = "n") |>
    transmute(metric = paste("Quadrant", .data$quadrant), .data$n,
              denominator = nrow(analysis))
) |>
  mutate(percent = 100 * .data$n / .data$denominator, cutoff = cutoff)
write_csv(summary_table, file.path(out_dir, "Supplementary_Figure_S1B_statistics.csv"))

long <- analysis |>
  select(.data$ModelID, .data$RIPK3, .data$NLRP3) |>
  pivot_longer(c("RIPK3", "NLRP3"), names_to = "gene", values_to = "expression") |>
  mutate(low = .data$expression < cutoff)
set.seed(20260814)
p_strip <- ggplot(long, aes(.data$expression, .data$gene, color = .data$low)) +
  geom_jitter(height = 0.16, width = 0, size = 1.2, alpha = 0.6) +
  geom_vline(xintercept = cutoff, linetype = "dashed", color = "#1F77B4") +
  scale_color_manual(values = c(`TRUE` = "#1F77B4", `FALSE` = "#D62728")) +
  labs(
    title = "RIPK3 and NLRP3 expression in human tumor cell-line models",
    subtitle = paste(metadata[["release"]], "; default expression entries; cutoff 0.5"),
    x = "Expression, log2(TPM+1)", y = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

p_scatter <- ggplot(analysis, aes(.data$RIPK3, .data$NLRP3)) +
  geom_point(alpha = 0.45, size = 1.5, color = "grey30") +
  geom_vline(xintercept = cutoff, linetype = "dashed", color = "#1F77B4") +
  geom_hline(yintercept = cutoff, linetype = "dashed", color = "#1F77B4") +
  labs(
    title = "RIPK3 versus NLRP3 expression",
    subtitle = sprintf("Human tumor cell-line models; n=%d", nrow(analysis)),
    x = "RIPK3, log2(TPM+1)", y = "NLRP3, log2(TPM+1)"
  ) +
  theme_classic(base_size = 14)

ggsave(file.path(out_dir, "Supplementary_Figure_S1B_strip.png"),
       p_strip, width = 8, height = 3.8, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Supplementary_Figure_S1B_scatter.png"),
       p_scatter, width = 6, height = 5.5, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Supplementary_Figure_S1B_scatter.tiff"),
       p_scatter, width = 6, height = 5.5, dpi = 600,
       compression = "lzw", bg = "white")

provenance <- list(
  depmap_release = metadata[["release"]],
  depmap_release_doi = metadata[["release_doi"]],
  download_date = metadata[["download_date"]],
  expression_sha256 = observed_expression_sha256,
  model_sha256 = observed_model_sha256,
  default_entry_field = default_flag,
  model_filters = list(
    ModelType = "Cell Line", TissueOrigin = "Human",
    OncotreePrimaryDisease = "non-empty"
  ),
  n_models = nrow(analysis),
  cutoff_log2_tpm_plus_1 = cutoff
)
write_json(
  provenance,
  file.path(out_dir, "Supplementary_Figure_S1B_provenance.json"),
  pretty = TRUE, auto_unbox = TRUE
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
