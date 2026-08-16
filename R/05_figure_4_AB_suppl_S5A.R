
# ============================================================
# FIGURE 4A-B, SUPPLEMENTARY FIGURE S5A, AND FIGURE 5A
#
# This script reproduces:
# - Figure 4A: UMAP of annotated T-cell clusters
# - Figure 4B: Cluster composition across WT, KO1, and KO2
# - Supplementary Figure S5A: UMAPs split by sample
# - Figure 5A: repeated CD3/CD28-stimulation analysis
#
# ------------------------------------------------------------
# DATA REQUIREMENTS (IMPORTANT)
#
# The version-controlled inputs are deterministic TSV.gz exports of the CD3+
# cell sheets from the four source workbooks. Their source and output SHA-256
# values are recorded in data/experimental/experimental_data_manifest.tsv.
#
# ------------------------------------------------------------
# OUTPUT
#
# The script generates figures and tables corresponding to:
# - Figure 4A-B
# - Supplementary Figure S5A
# - Figure 5A
#
# Output figures are saved to:
#
#   figures/
#
# ------------------------------------------------------------
# NOTE
#
# This script processes BD Rhapsody targeted T-cell panel data,
# including WT, KO1, and KO2 samples, and performs:
#
# - quality control filtering
# - normalization and scaling
# - dimensionality reduction (PCA, UMAP)
# - clustering and manual annotation
# - cluster composition analysis
# - marker gene identification
#
# C10 is retained in the complete object and marker/QC exports. It is excluded post-clustering from the
# C0-C9 state summaries and transferred signatures used in the manuscript.
# This is an explicit analysis decision, not a sample-level QC exclusion. The
# repeated-stimulation dataset is clustered and annotated independently; its
# C0-C5 labels are not mapped onto the tumor-co-culture C0-C5 labels.
#
# The script generates the following figures:
# - Figure 4A-B (UMAP and cluster composition)
# - Supplementary Figure S5A (UMAP per sample)
# - Figure 5A (repeated-stimulation analysis)
#
# A separate section performs the repeated-stimulation analysis using the
# TCR CD3+ subset (5,662 cells); the three CD3-negative cells in the source
# workbook's "all cells" sheet are not included.

# -----------------------------
# Install/load packages
# -----------------------------
install_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) return(invisible(TRUE))
  stop("Missing required packages: ", paste(missing, collapse = ", "),
       ". Install and record the required versions before running.")
}

needed <- c(
  "Matrix", "Seurat", "dplyr", "tidyr", "tibble",
  "ggplot2", "patchwork", "data.table", "R.utils", "scales", "openxlsx"
)
install_if_missing(needed)

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
  library(openxlsx)
  library(grid)
})

# -----------------------------
# -----------------------------
# Input files
# -----------------------------
FILES <- c(
  WT  = file.path("data", "experimental", "singlecell", "WT_targeted_counts.tsv.gz"),
  KO1 = file.path("data", "experimental", "singlecell", "KO1_targeted_counts.tsv.gz"),
  KO2 = file.path("data", "experimental", "singlecell", "KO2_targeted_counts.tsv.gz")
)

# -----------------------------
# Output directory. All outputs from this script are staged and are exposed
# only after both marker/annotation guards have passed.
# -----------------------------
FINAL_FIG_DIR <- "figures"
DIAGNOSTIC_DIR <- file.path("results", "targeted-singlecell-diagnostics")
FIG_DIR <- tempfile(pattern = "R05-staging-", tmpdir = tempdir())
if (!dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)) {
  stop("Could not create the R/05 staging directory: ", FIG_DIR)
}

promote_staged_outputs <- function(staging_dir, final_dir) {
  staged_paths <- list.files(
    staging_dir,
    all.files = TRUE,
    full.names = TRUE,
    recursive = TRUE,
    include.dirs = FALSE,
    no.. = TRUE
  )
  if (!length(staged_paths)) stop("R/05 produced no staged outputs.")

  relative_paths <- substring(staged_paths, nchar(staging_dir) + 2L)
  for (index in seq_along(staged_paths)) {
    target_path <- file.path(final_dir, relative_paths[[index]])
    dir.create(dirname(target_path), recursive = TRUE, showWarnings = FALSE)
    incoming_path <- tempfile(
      pattern = paste0(".", basename(target_path), ".incoming-"),
      tmpdir = dirname(target_path)
    )
    copied <- file.copy(
      staged_paths[[index]], incoming_path,
      overwrite = FALSE, copy.mode = TRUE, copy.date = TRUE
    )
    if (!copied || !identical(
      unname(tools::md5sum(staged_paths[[index]])),
      unname(tools::md5sum(incoming_path))
    )) {
      unlink(incoming_path, force = TRUE)
      stop("Could not stage a verified copy for publication: ", target_path)
    }

    # Atomic replacement on POSIX. file.rename() cannot replace an existing
    # file on some Windows installations, so use a guarded copy fallback there.
    if (!file.rename(incoming_path, target_path)) {
      copied <- file.copy(
        incoming_path, target_path,
        overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE
      )
      unlink(incoming_path, force = TRUE)
      if (!copied) stop("Could not publish validated output: ", target_path)
    }
  }
  invisible(relative_paths)
}

QC_Q_FEATURE_LOW  <- 0.02
QC_Q_FEATURE_HIGH <- 0.98
QC_Q_COUNT_HIGH   <- 0.99

NORM_SCALE_FACTOR <- 1e4
REGRESS_VARS      <- c("nCount_RNA")
VARSEL_METHOD     <- "vst"
VARSEL_NFEATURES  <- function(n_genes) min(2000, n_genes)

PCA_NPCS          <- 30
NEIGHBOR_DIMS     <- 1:20
FINAL_RES         <- 0.6
UMAP_DIMS         <- 1:20
UMAP_N_NEIGHBORS  <- 30
UMAP_MIN_DIST     <- 0.3

# -----------------------------
# Cluster labels
# -----------------------------
cluster_labels_short <- c(
  "0"  = "C0",
  "1"  = "C1",
  "2"  = "C2",
  "3"  = "C3",
  "4"  = "C4",
  "5"  = "C5",
  "6"  = "C6",
  "7"  = "C7",
  "8"  = "C8",
  "9"  = "C9",
  "10" = "C10"
)

cluster_labels_full <- c(
  "C0"  = "C0 KLRB1+ innate-like activated",
  "C1"  = "C1 activated CD4 checkpoint-high",
  "C2"  = "C2 CXCR5+ type 2 helper-like",
  "C3"  = "C3 CXCR6+ cytotoxic activated",
  "C4"  = "C4 ZNF683+ cytotoxic CD8",
  "C5"  = "C5 γδ T",
  "C6"  = "C6 CXCL13-associated cycling T-cell state",
  "C7"  = "C7 cycling effector-like",
  "C8"  = "C8 IL9hi TH9-like",
  "C9"  = "C9 TCF7hi stem-like/early-memory",
  "C10" = "C10 post-clustering excluded cluster"
)

# -----------------------------
# Colors
# -----------------------------
cols_samples <- c(
  WT  = "#0000FF",
  KO1 = "#33CC00",
  KO2 = "#FB0207"
)

cluster_cols <- c(
  "C0 KLRB1+ innate-like activated"         = "#F08A80",
  "C1 activated CD4 checkpoint-high"        = "#D98C00",
  "C2 CXCR5+ type 2 helper-like"            = "#B8A500",
  "C3 CXCR6+ cytotoxic activated"           = "#6CB400",
  "C4 ZNF683+ cytotoxic CD8"                = "#20B95A",
  "C5 γδ T"                                 = "#16B8A8",
  "C6 CXCL13-associated cycling T-cell state" = "#29AFC4",
  "C7 cycling effector-like"                = "#2A96E6",
  "C8 IL9hi TH9-like"                       = "#A27AE8",
  "C9 TCF7hi stem-like/early-memory"        = "#D865D8",
  "C10 post-clustering excluded cluster"    = "#E75AA2"
)

# -----------------------------
# Utilities
# -----------------------------
save_plot <- function(p, fn, w = 7, h = 6, dpi = 600) {
  png_path <- file.path(FIG_DIR, fn)
  ggsave(
    filename = png_path,
    plot = p,
    width = w,
    height = h,
    dpi = dpi,
    bg = "white"
  )
  tiff_path <- file.path(
    FIG_DIR,
    paste0(tools::file_path_sans_ext(basename(fn)), ".tiff")
  )
  ggsave(
    filename = tiff_path,
    plot = p,
    width = w,
    height = h,
    dpi = dpi,
    compression = "lzw",
    bg = "white"
  )
}

read_targeted_count_matrix <- function(path) {
  df <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  if (ncol(df) < 2L || names(df)[1] != "cell_index") {
    stop("Unexpected targeted-count table schema: ", path)
  }

  cell_ids <- trimws(as.character(df[[1]]))
  if (any(is.na(cell_ids) | cell_ids == "") || anyDuplicated(cell_ids)) {
    stop("Missing or duplicate cell identifier(s) in: ", path)
  }
  expr_df <- df[, -1, drop = FALSE]
  gene_names <- trimws(colnames(expr_df))
  if (any(is.na(gene_names) | gene_names == "") || anyDuplicated(gene_names)) {
    stop("Missing or duplicate gene name(s) in: ", path)
  }
  expr_mat_cell_gene <- as.matrix(expr_df)
  storage.mode(expr_mat_cell_gene) <- "numeric"
  if (any(!is.finite(expr_mat_cell_gene))) {
    stop("The count matrix contains missing or non-numeric value(s): ", path)
  }
  if (any(expr_mat_cell_gene < 0) || any(expr_mat_cell_gene %% 1 != 0)) {
    stop("The count matrix must contain non-negative integer counts: ", path)
  }
  rownames(expr_mat_cell_gene) <- cell_ids

  counts <- t(expr_mat_cell_gene)
  rownames(counts) <- gene_names
  colnames(counts) <- cell_ids

  counts <- Matrix::Matrix(counts, sparse = TRUE)

  list(
    counts = counts,
    source_subset = "CD3+ cells",
    n_genes = nrow(counts),
    n_cells = ncol(counts)
  )
}

# -----------------------------
# Read files and create objects
# -----------------------------
set.seed(12345)
objs <- list()

for (nm in names(FILES)) {
  fp <- FILES[[nm]]
  if (!file.exists(fp)) stop("Missing input file: ", fp)

  dat <- read_targeted_count_matrix(fp)
  counts <- dat$counts

  message("Loaded ", nm,
          " | source subset=", dat$source_subset,
          " | genes=", dat$n_genes,
          " | cells=", dat$n_cells)

  obj_i <- CreateSeuratObject(
    counts = counts,
    project = nm,
    min.cells = 0,
    min.features = 0
  )

  obj_i$sample <- nm
  obj_i$condition <- ifelse(nm == "WT", "WT", "KO")
  objs[[nm]] <- obj_i
}

# -----------------------------
# Merge + Seurat v5-safe layers
# -----------------------------
obj <- merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs), project = "WT_vs_KO")
obj <- JoinLayers(obj, assay = "RNA")

if (!"sample" %in% colnames(obj@meta.data)) {
  if ("orig.ident" %in% colnames(obj@meta.data)) {
    obj$sample <- as.character(obj$orig.ident)
  } else {
    obj$sample <- sub("_.*$", "", colnames(obj))
  }
}
obj$sample <- factor(as.character(obj$sample), levels = c("WT", "KO1", "KO2"))
obj$condition <- ifelse(obj$sample == "WT", "WT", "KO")

# -----------------------------
# QC
# -----------------------------
rna_counts <- LayerData(obj, assay = "RNA", layer = "counts")
obj$nCount_RNA <- Matrix::colSums(rna_counts)
obj$nFeature_RNA <- Matrix::colSums(rna_counts > 0)

meta_before <- obj@meta.data %>%
  tibble::rownames_to_column("cell")

if (!"sample" %in% colnames(meta_before)) {
  if ("orig.ident" %in% colnames(meta_before)) {
    meta_before$sample <- as.character(meta_before$orig.ident)
  } else {
    meta_before$sample <- sub("_.*$", "", meta_before$cell)
  }
}
meta_before$sample <- factor(as.character(meta_before$sample), levels = c("WT", "KO1", "KO2"))

qc_thresholds <- meta_before %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(
    nFeature_low  = quantile(nFeature_RNA, QC_Q_FEATURE_LOW,  na.rm = TRUE),
    nFeature_high = quantile(nFeature_RNA, QC_Q_FEATURE_HIGH, na.rm = TRUE),
    nCount_high   = quantile(nCount_RNA,   QC_Q_COUNT_HIGH,   na.rm = TRUE),
    .groups = "drop"
  )

meta_qc <- meta_before %>%
  dplyr::left_join(qc_thresholds, by = "sample") %>%
  dplyr::mutate(
    pass_qc = (nFeature_RNA >= nFeature_low) &
      (nFeature_RNA <= nFeature_high) &
      (nCount_RNA   <= nCount_high)
  )

qc_summary <- meta_qc %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(
    input_cells = dplyr::n(),
    retained_cells = sum(pass_qc),
    excluded_cells = sum(!pass_qc),
    nFeature_low = dplyr::first(nFeature_low),
    nFeature_high = dplyr::first(nFeature_high),
    nCount_high = dplyr::first(nCount_high),
    .groups = "drop"
  )
data.table::fwrite(
  qc_summary,
  file.path(FIG_DIR, "Supplementary_Table_S5_cell_filtering_QC.tsv"),
  sep = "\t"
)

cells_keep <- meta_qc$cell[meta_qc$pass_qc]
obj <- subset(obj, cells = cells_keep)

if (!"sample" %in% colnames(obj@meta.data)) {
  if ("orig.ident" %in% colnames(obj@meta.data)) {
    obj$sample <- as.character(obj$orig.ident)
  } else {
    obj$sample <- sub("_.*$", "", colnames(obj))
  }
}
obj$sample <- factor(as.character(obj$sample), levels = c("WT", "KO1", "KO2"))
obj$condition <- ifelse(obj$sample == "WT", "WT", "KO")

rna_counts <- LayerData(obj, assay = "RNA", layer = "counts")
obj$nCount_RNA <- Matrix::colSums(rna_counts)
obj$nFeature_RNA <- Matrix::colSums(rna_counts > 0)

# -----------------------------
# Normalize / variable features / scale
# -----------------------------
obj <- NormalizeData(
  obj,
  normalization.method = "LogNormalize",
  scale.factor = NORM_SCALE_FACTOR,
  verbose = FALSE
)

nfeatures_use <- if (is.function(VARSEL_NFEATURES)) {
  VARSEL_NFEATURES(nrow(obj))
} else {
  VARSEL_NFEATURES
}
nfeatures_use <- as.integer(nfeatures_use[1])
nfeatures_use <- min(nfeatures_use, nrow(obj))
nfeatures_use <- max(nfeatures_use, 50L)

obj <- FindVariableFeatures(
  obj,
  selection.method = VARSEL_METHOD,
  nfeatures = nfeatures_use,
  verbose = FALSE
)

vars_regress_use <- REGRESS_VARS[REGRESS_VARS %in% colnames(obj@meta.data)]
obj <- ScaleData(
  obj,
  vars.to.regress = vars_regress_use,
  features = rownames(obj),
  verbose = FALSE
)

# -----------------------------
# PCA / clustering / UMAP
# -----------------------------
pca_features <- VariableFeatures(obj)
if (length(pca_features) < 2) stop("Too few variable features for PCA.")

npcs_use <- min(as.integer(PCA_NPCS), length(pca_features))
npcs_use <- max(npcs_use, 2L)

obj <- RunPCA(
  obj,
  npcs = npcs_use,
  features = pca_features,
  verbose = FALSE
)

dims_nn <- NEIGHBOR_DIMS[NEIGHBOR_DIMS <= npcs_use]
dims_umap <- UMAP_DIMS[UMAP_DIMS <= npcs_use]

if (length(dims_nn) < 2) stop("Too few PCA dimensions for FindNeighbors.")
if (length(dims_umap) < 2) stop("Too few PCA dimensions for RunUMAP.")

obj <- FindNeighbors(
  obj,
  dims = dims_nn,
  verbose = FALSE
)

obj <- FindClusters(
  obj,
  resolution = FINAL_RES,
  algorithm = 1,
  random.seed = 12345,
  verbose = FALSE
)

obj$seurat_clusters <- as.character(Idents(obj))
observed_cluster_ids <- sort(unique(obj$seurat_clusters))
expected_cluster_ids <- names(cluster_labels_short)
if (!setequal(observed_cluster_ids, expected_cluster_ids)) {
  stop(
    "The frozen C0-C10 annotation requires clusters 0-10; observed: ",
    paste(observed_cluster_ids, collapse = ", ")
  )
}

obj <- RunUMAP(
  obj,
  dims = dims_umap,
  n.neighbors = UMAP_N_NEIGHBORS,
  min.dist = UMAP_MIN_DIST,
  umap.method = "uwot",
  seed.use = 12345,
  verbose = FALSE
)

# -----------------------------
# Cluster annotations for final figures
# -----------------------------
obj$cluster_short <- factor(
  unname(cluster_labels_short[as.character(obj$seurat_clusters)]),
  levels = unname(cluster_labels_short)
)

obj$cluster_annot <- factor(
  unname(cluster_labels_full[as.character(obj$cluster_short)]),
  levels = unname(cluster_labels_full[paste0("C", 0:10)])
)

# Keep C10 in the complete object and QC outputs, but exclude it from the
# manuscript's C0-C9 state summaries and transferred signatures. The complete
# counts below make this decision auditable.
obj_biological <- subset(obj, subset = cluster_short != "C10")

if (ncol(obj_biological) + sum(obj$cluster_short == "C10") != ncol(obj)) {
  stop("C10 exclusion check failed.")
}

message(
  "Cells retained for biological analyses: ", ncol(obj_biological),
  "; C10 cells excluded: ", sum(obj$cluster_short == "C10")
)

# -----------------------------
# FIGURE 1: UMAP by sample
# -----------------------------
p_umap_sample <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "sample",
  cols = cols_samples,
  pt.size = 0.4
) +
  ggtitle("UMAP by sample") +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 22),
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 14, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 15)
  )

save_plot(p_umap_sample, "Supplementary_Figure_S5A_UMAP_by_sample.png", w = 9, h = 7)

# -----------------------------
# FIGURE 2: 4-panel UMAP
# -----------------------------
p_WT <- DimPlot(
  subset(obj, subset = sample == "WT"),
  reduction = "umap",
  group.by = "sample",
  cols = cols_samples,
  pt.size = 0.4
) +
  ggtitle("WT") +
  NoLegend() +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11, color = "black")
  )

p_KO1 <- DimPlot(
  subset(obj, subset = sample == "KO1"),
  reduction = "umap",
  group.by = "sample",
  cols = cols_samples,
  pt.size = 0.4
) +
  ggtitle("KO1") +
  NoLegend() +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11, color = "black")
  )

p_KO2 <- DimPlot(
  subset(obj, subset = sample == "KO2"),
  reduction = "umap",
  group.by = "sample",
  cols = cols_samples,
  pt.size = 0.4
) +
  ggtitle("KO2") +
  NoLegend() +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11, color = "black")
  )

p_all <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "sample",
  cols = cols_samples,
  pt.size = 0.4
) +
  ggtitle("All samples") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

p_4panel <- p_WT | p_KO1 | p_KO2 | p_all

save_plot(p_4panel, "Figure_4A_UMAP_4panel.png", w = 17, h = 5)

# -----------------------------
# FIGURE 3: UMAP clusters clean
# -----------------------------
umap_coordinates <- as.data.frame(Embeddings(obj_biological, "umap"))
colnames(umap_coordinates)[1:2] <- c("umap_1", "umap_2")
centers <- obj_biological@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::bind_cols(umap_coordinates) %>%
  dplyr::group_by(cluster_short) %>%
  dplyr::summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2),
    .groups = "drop"
  )

p_clusters_clean <- DimPlot(
  obj_biological,
  reduction = "umap",
  group.by = "cluster_annot",
  cols = cluster_cols,
  label = FALSE,
  pt.size = 0.4
) +
  geom_text(
    data = centers,
    aes(x = umap_1, y = umap_2, label = cluster_short),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 6
  ) +
  ggtitle("UMAP by annotated clusters") +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 22),
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 14, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

save_plot(p_clusters_clean, "Figure_4A_UMAP_clusters_annotated.png", w = 15, h = 10)

# -----------------------------
# FIGURE 4: Classic stacked barplot
# -----------------------------
md_all <- obj@meta.data %>%
  tibble::rownames_to_column("cell")

if (!"sample" %in% colnames(md_all)) {
  if ("orig.ident" %in% colnames(md_all)) {
    md_all$sample <- as.character(md_all$orig.ident)
  } else {
    md_all$sample <- sub("_.*$", "", md_all$cell)
  }
}
md_all$sample <- factor(md_all$sample, levels = c("WT", "KO1", "KO2"))

md_all$cluster_short <- factor(
  unname(cluster_labels_short[as.character(md_all$seurat_clusters)]),
  levels = paste0("C", 0:10)
)

md_all$cluster_full <- factor(
  unname(cluster_labels_full[as.character(md_all$cluster_short)]),
  levels = unname(cluster_labels_full[paste0("C", 0:10)])
)

md <- md_all %>%
  dplyr::filter(cluster_short != "C10") %>%
  dplyr::mutate(
    cluster_short = droplevels(cluster_short),
    cluster_full = droplevels(cluster_full)
  )

df_stack <- md %>%
  dplyr::group_by(sample, cluster_full) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(percent = 100 * n / sum(n)) %>%
  dplyr::ungroup()

p_bar_classic <- ggplot(df_stack, aes(x = sample, y = percent, fill = cluster_full)) +
  geom_col(color = "white", linewidth = 1.2, width = 0.88) +
  scale_fill_manual(values = cluster_cols) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  labs(
    title = "Descriptive cluster composition by sample",
    x = "Sample",
    y = "% of C0-C9 QC-passing cells",
    fill = NULL
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 24),
    axis.title = element_text(face = "bold", size = 20),
    axis.text = element_text(face = "bold", size = 16, color = "black"),
    legend.text = element_text(size = 15),
    legend.key.size = unit(0.7, "cm"),
    legend.position = "right"
  )

save_plot(p_bar_classic, "Figure_4B_cluster_composition_optional.png", w = 14, h = 10)

# -----------------------------
# TABLE 1: Cell counts by cluster and sample
# -----------------------------
# QC table retaining C10 for complete reporting.
cluster_cell_counts_all_qc <- md_all %>%
  dplyr::group_by(cluster_short, sample) %>%
  dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = sample,
    values_from = n_cells,
    values_fill = 0
  ) %>%
  dplyr::arrange(cluster_short) %>%
  dplyr::mutate(Total = WT + KO1 + KO2)

write.csv(
  cluster_cell_counts_all_qc,
  file.path(FIG_DIR, "Supplementary_Table_S5_cell_counts_all_QC_including_C10.csv"),
  row.names = FALSE
)

# Biological table used for Figure 4B and downstream signatures (C10 excluded).
cluster_cell_counts <- md %>%
  dplyr::group_by(cluster_short, sample) %>%
  dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = sample,
    values_from = n_cells,
    values_fill = 0
  ) %>%
  dplyr::arrange(cluster_short) %>%
  dplyr::mutate(Total = WT + KO1 + KO2)

write.csv(
  cluster_cell_counts,
  file.path(FIG_DIR, "Supplementary_Table_S5_cell_counts_by_cluster.csv"),
  row.names = FALSE
)

# -----------------------------
# Exploratory cell-level marker rankings. The C0-C9-only run reproduces the
# comparison population used to freeze the transferred signatures. A separate
# all-cluster run retains C10-versus-rest evidence for reviewing its
# post-clustering exclusion. These tests do not constitute biological-replicate
# inference between experimental groups.
# -----------------------------
DefaultAssay(obj) <- "RNA"
DefaultAssay(obj_biological) <- "RNA"

markers_all <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
if (nrow(markers_all) == 0L) {
  stop("All-cluster marker table is empty; C10 exclusion cannot be reviewed.")
}
marker_counts_all <- markers_all %>%
  dplyr::mutate(cluster = as.character(.data$cluster)) %>%
  dplyr::count(.data$cluster, name = "n_markers")
expected_all_clusters <- as.character(0:10)
if (!setequal(marker_counts_all$cluster, expected_all_clusters) ||
    any(marker_counts_all$n_markers[
      match(expected_all_clusters, marker_counts_all$cluster)
    ] < 20L)) {
  stop("Every cluster 0-10 must have at least 20 ranked markers for review.")
}
data.table::fwrite(
  markers_all,
  file.path(FIG_DIR, "Supplementary_Table_S5_all_cluster_markers_including_C10.tsv.gz"),
  sep = "\t"
)
data.table::fwrite(
  markers_all %>% dplyr::filter(as.character(.data$cluster) == "10"),
  file.path(FIG_DIR, "Supplementary_Table_S5_C10_markers.tsv"),
  sep = "\t"
)

markers <- FindAllMarkers(
  obj_biological,
  only.pos = TRUE,
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

if (nrow(markers) > 0) {
  marker_counts <- markers %>%
    dplyr::mutate(cluster = as.character(.data$cluster)) %>%
    dplyr::count(.data$cluster, name = "n_markers")
  expected_marker_clusters <- as.character(0:9)
  if (!setequal(marker_counts$cluster, expected_marker_clusters) ||
      any(marker_counts$n_markers[match(expected_marker_clusters, marker_counts$cluster)] < 20L)) {
    stop("Every C0-C9 cluster must have at least 20 ranked markers for review.")
  }
  fc_col <- if ("avg_log2FC" %in% colnames(markers)) {
    "avg_log2FC"
  } else if ("avg_logFC" %in% colnames(markers)) {
    "avg_logFC"
  } else {
    NULL
  }

  if (!is.null(fc_col)) {
    top_markers <- markers %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_max(order_by = .data[[fc_col]], n = 20, with_ties = FALSE) %>%
      dplyr::ungroup()

    frozen_signatures <- data.table::fread(
      file.path("resources", "CAR_T_state_signatures.csv"), data.table = FALSE
    ) %>%
      dplyr::transmute(
        cluster = paste0("C", as.character(.data$cluster)),
        gene = toupper(trimws(as.character(.data$gene)))
      )
    current_membership <- top_markers %>%
      dplyr::transmute(
        cluster = paste0("C", as.character(.data$cluster)),
        gene = toupper(trimws(as.character(.data$gene)))
      )
    signature_concordance <- dplyr::bind_rows(lapply(paste0("C", 0:9), function(cl) {
      frozen_genes <- frozen_signatures$gene[frozen_signatures$cluster == cl]
      current_genes <- current_membership$gene[current_membership$cluster == cl]
      dplyr::tibble(
        cluster = cl,
        n_frozen = length(frozen_genes),
        n_current = length(current_genes),
        n_overlap = length(intersect(frozen_genes, current_genes)),
        exact_membership_match = setequal(frozen_genes, current_genes)
      )
    }))
    data.table::fwrite(
      signature_concordance,
      file.path(FIG_DIR, "Supplementary_Table_S5_signature_concordance.tsv"),
      sep = "\t"
    )
    if (any(!signature_concordance$exact_membership_match)) {
      mismatched_clusters <- signature_concordance$cluster[
        !signature_concordance$exact_membership_match
      ]
      signature_differences <- dplyr::bind_rows(lapply(
        mismatched_clusters,
        function(cl) {
          frozen_genes <- frozen_signatures$gene[frozen_signatures$cluster == cl]
          current_genes <- current_membership$gene[current_membership$cluster == cl]
          frozen_only <- sort(setdiff(frozen_genes, current_genes))
          current_only <- sort(setdiff(current_genes, frozen_genes))
          dplyr::bind_rows(
            dplyr::tibble(
              cluster = rep(cl, length(frozen_only)),
              membership = rep("frozen_only", length(frozen_only)),
              gene = frozen_only
            ),
            dplyr::tibble(
              cluster = rep(cl, length(current_only)),
              membership = rep("current_only", length(current_only)),
              gene = current_only
            )
          )
        }
      ))
      dir.create(DIAGNOSTIC_DIR, recursive = TRUE, showWarnings = FALSE)
      data.table::fwrite(
        signature_concordance,
        file.path(DIAGNOSTIC_DIR, "Supplementary_Table_S5_signature_concordance.tsv"),
        sep = "\t"
      )
      data.table::fwrite(
        top_markers,
        file.path(DIAGNOSTIC_DIR, "Supplementary_Table_S5_current_top20_markers.tsv"),
        sep = "\t"
      )
      data.table::fwrite(
        frozen_signatures,
        file.path(DIAGNOSTIC_DIR, "Supplementary_Table_S5_frozen_signatures.tsv"),
        sep = "\t"
      )
      data.table::fwrite(
        signature_differences,
        file.path(DIAGNOSTIC_DIR, "Supplementary_Table_S5_signature_membership_differences.tsv"),
        sep = "\t"
      )
      message("C0-C9 frozen-signature concordance:")
      print(as.data.frame(signature_concordance), row.names = FALSE)
      for (cl in mismatched_clusters) {
        cluster_differences <- signature_differences[
          signature_differences$cluster == cl,
          ,
          drop = FALSE
        ]
        for (membership in c("frozen_only", "current_only")) {
          genes <- cluster_differences$gene[
            cluster_differences$membership == membership
          ]
          message(
            cl, " ", membership, ": ",
            if (length(genes)) paste(genes, collapse = ", ") else "<none>"
          )
        }
      }
      stop(
        "Current C0-C9 top-20 markers do not match the frozen signatures; ",
        "manual marker/label review is required before figure release. ",
        "Diagnostic tables were written to ", DIAGNOSTIC_DIR, "."
      )
    }

    wb <- openxlsx::createWorkbook()

    # summary sheet
    openxlsx::addWorksheet(wb, "All_clusters_top20")
    openxlsx::writeData(wb, "All_clusters_top20", top_markers)

    clusters <- sort(unique(top_markers$cluster))

    for (cl in clusters) {
      df_cl <- top_markers %>% dplyr::filter(cluster == cl)
      sheet_name <- paste0("Cluster_", cl)
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeData(wb, sheet_name, df_cl)
    }

    openxlsx::saveWorkbook(
      wb,
      file.path(FIG_DIR, "Supplementary_Table_S5_top_markers_per_cluster.xlsx"),
      overwrite = TRUE
    )
  } else {
    stop("Could not find avg_log2FC or avg_logFC column in markers table.")
  }
} else {
  stop("Marker table is empty; cluster annotations cannot be release-validated.")
}

# -----------------------------

# -----------------------------
# FIGURE 5: Faceted barplot
# -----------------------------
df_bar <- md %>%
  dplyr::group_by(sample, cluster_short) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(frac = n / sum(n)) %>%
  dplyr::ungroup()

p_bar_facet <- ggplot(df_bar, aes(x = cluster_short, y = frac, fill = sample)) +
  geom_col(width = 0.8) +
  facet_grid(sample ~ ., switch = "y") +
  scale_fill_manual(values = cols_samples) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    position = "right"
  ) +
  labs(
    x = "Cluster",
    y = "Fraction of cells"
  ) +
  theme_bw(base_size = 16) +
  theme(
    strip.text.y.left = element_text(face = "bold", size = 16, angle = 90),
    strip.background = element_rect(fill = "grey85", color = "grey30", linewidth = 1),
    strip.placement = "outside",

    axis.title = element_text(face = "bold", size = 16),
    axis.title.y.right = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 0.5),  #

    axis.text.x = element_text(
      angle = 0,
      vjust = 0.5,
      hjust = 0.5,   #
      size = 12,
      color = "black"
    ),
    axis.text.y = element_text(size = 12, color = "black"),

    legend.position = "none",
    panel.spacing = unit(0.15, "lines")
  )

save_plot(p_bar_facet, "Figure_4B_cluster_fraction_facet.png", w = 4.8, h = 10)













# ============================================================
# Repeated-stimulation analysis of BD Rhapsody targeted T-cell panel data.
# TCR cluster numbers are dataset-specific and are not mapped onto the
# independently clustered CAR-T dataset.
# ============================================================

# -----------------------------
# Install/load packages
# -----------------------------
install_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) return(invisible(TRUE))
  stop("Missing required packages: ", paste(missing, collapse = ", "),
       ". Install and record the required versions before running.")
}

needed <- c(
  "Matrix", "Seurat", "dplyr", "tidyr", "tibble",
  "ggplot2", "patchwork", "data.table", "R.utils", "scales", "openxlsx"
)
install_if_missing(needed)

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
  library(openxlsx)
  library(grid)
})

# -----------------------------
# User parameters
# -----------------------------
# -----------------------------
# Input file
# -----------------------------
FILE_TCR <- file.path(
  "data", "experimental", "singlecell", "TCR_targeted_counts.tsv.gz"
)
SOURCE_SUBSET <- "CD3+ cells"

# -----------------------------
# Continue using the run-level staging directory initialized above.
# -----------------------------

QC_Q_FEATURE_LOW  <- 0.02
QC_Q_FEATURE_HIGH <- 0.98
QC_Q_COUNT_HIGH   <- 0.99

NORM_SCALE_FACTOR <- 1e4
REGRESS_VARS      <- c("nCount_RNA")
VARSEL_METHOD     <- "vst"
VARSEL_NFEATURES  <- function(n_genes) min(2000, n_genes)

PCA_NPCS         <- 30
NEIGHBOR_DIMS    <- 1:20
FINAL_RES        <- 0.6
UMAP_DIMS        <- 1:20
UMAP_N_NEIGHBORS <- 30
UMAP_MIN_DIST    <- 0.3

# -----------------------------
# Colors
# -----------------------------
cluster_palette <- c(
  "#F08A80", "#D98C00", "#B8A500", "#6CB400", "#20B95A",
  "#16B8A8", "#29AFC4", "#2A96E6", "#A27AE8", "#D865D8",
  "#E75AA2", "#9E9E9E", "#8D6E63", "#5C6BC0", "#26A69A",
  "#EC407A", "#7CB342", "#FF7043", "#AB47BC", "#78909C"
)

# -----------------------------
# Utilities
# -----------------------------
save_plot <- function(p, fn, w = 7, h = 6, dpi = 600) {
  png_path <- file.path(FIG_DIR, fn)
  ggsave(
    filename = png_path,
    plot = p,
    width = w,
    height = h,
    dpi = dpi,
    bg = "white"
  )
  tiff_path <- file.path(
    FIG_DIR,
    paste0(tools::file_path_sans_ext(basename(fn)), ".tiff")
  )
  ggsave(
    filename = tiff_path,
    plot = p,
    width = w,
    height = h,
    dpi = dpi,
    compression = "lzw",
    bg = "white"
  )
}

read_targeted_count_matrix <- function(path) {
  df <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  if (ncol(df) < 2L || names(df)[1] != "cell_index") {
    stop("Unexpected targeted-count table schema: ", path)
  }
  cell_ids <- trimws(as.character(df[[1]]))
  if (any(is.na(cell_ids) | cell_ids == "") || anyDuplicated(cell_ids)) {
    stop("Missing or duplicate cell identifier(s) in: ", path)
  }
  expr_df <- df[, -1, drop = FALSE]
  gene_names <- trimws(colnames(expr_df))
  if (any(is.na(gene_names) | gene_names == "") || anyDuplicated(gene_names)) {
    stop("Missing or duplicate gene name(s) in: ", path)
  }
  expr_mat_cell_gene <- as.matrix(expr_df)
  storage.mode(expr_mat_cell_gene) <- "numeric"
  if (any(!is.finite(expr_mat_cell_gene))) {
    stop("The count matrix contains missing or non-numeric value(s): ", path)
  }
  if (any(expr_mat_cell_gene < 0) || any(expr_mat_cell_gene %% 1 != 0)) {
    stop("The count matrix must contain non-negative integer counts: ", path)
  }
  rownames(expr_mat_cell_gene) <- cell_ids

  counts <- t(expr_mat_cell_gene)
  rownames(counts) <- gene_names
  colnames(counts) <- cell_ids
  counts <- Matrix::Matrix(counts, sparse = TRUE)

  list(
    counts = counts,
    source_subset = SOURCE_SUBSET,
    n_genes = nrow(counts),
    n_cells = ncol(counts)
  )
}

# -----------------------------
# Read TCR file
# -----------------------------
set.seed(12345)
if (!file.exists(FILE_TCR)) stop("Missing input file: ", FILE_TCR)

dat <- read_targeted_count_matrix(FILE_TCR)
counts <- dat$counts

message(
  "Loaded TCR | source subset=", dat$source_subset,
  " | genes=", dat$n_genes,
  " | cells=", dat$n_cells
)

obj <- CreateSeuratObject(
  counts = counts,
  project = "TCR",
  min.cells = 0,
  min.features = 0
)

obj$sample <- "TCR"
obj$source_subset <- dat$source_subset

obj <- JoinLayers(obj, assay = "RNA")

# -----------------------------
# QC
# -----------------------------
rna_counts <- LayerData(obj, assay = "RNA", layer = "counts")
obj$nCount_RNA   <- Matrix::colSums(rna_counts)
obj$nFeature_RNA <- Matrix::colSums(rna_counts > 0)

meta_before <- obj@meta.data %>%
  tibble::rownames_to_column("cell")

nFeature_low  <- quantile(meta_before$nFeature_RNA, QC_Q_FEATURE_LOW,  na.rm = TRUE)
nFeature_high <- quantile(meta_before$nFeature_RNA, QC_Q_FEATURE_HIGH, na.rm = TRUE)
nCount_high   <- quantile(meta_before$nCount_RNA,   QC_Q_COUNT_HIGH,   na.rm = TRUE)

meta_qc <- meta_before %>%
  dplyr::mutate(
    pass_qc = (nFeature_RNA >= nFeature_low) &
      (nFeature_RNA <= nFeature_high) &
      (nCount_RNA   <= nCount_high)
  )

tcr_qc_summary <- data.frame(
  sample = "TCR",
  input_cells = nrow(meta_qc),
  retained_cells = sum(meta_qc$pass_qc),
  excluded_cells = sum(!meta_qc$pass_qc),
  nFeature_low = nFeature_low,
  nFeature_high = nFeature_high,
  nCount_high = nCount_high
)
data.table::fwrite(
  tcr_qc_summary,
  file.path(FIG_DIR, "Figure_5A_TCR_cell_filtering_QC.tsv"),
  sep = "\t"
)

cells_keep <- meta_qc$cell[meta_qc$pass_qc]
obj <- subset(obj, cells = cells_keep)

rna_counts <- LayerData(obj, assay = "RNA", layer = "counts")
obj$nCount_RNA   <- Matrix::colSums(rna_counts)
obj$nFeature_RNA <- Matrix::colSums(rna_counts > 0)

# -----------------------------
# Normalize / variable features / scale
# -----------------------------
obj <- NormalizeData(
  obj,
  normalization.method = "LogNormalize",
  scale.factor = NORM_SCALE_FACTOR,
  verbose = FALSE
)

nfeatures_use <- if (is.function(VARSEL_NFEATURES)) {
  VARSEL_NFEATURES(nrow(obj))
} else {
  VARSEL_NFEATURES
}
nfeatures_use <- as.integer(nfeatures_use[1])
nfeatures_use <- min(nfeatures_use, nrow(obj))
nfeatures_use <- max(nfeatures_use, 50L)

obj <- FindVariableFeatures(
  obj,
  selection.method = VARSEL_METHOD,
  nfeatures = nfeatures_use,
  verbose = FALSE
)

vars_regress_use <- REGRESS_VARS[REGRESS_VARS %in% colnames(obj@meta.data)]
obj <- ScaleData(
  obj,
  vars.to.regress = vars_regress_use,
  features = rownames(obj),
  verbose = FALSE
)

# -----------------------------
# PCA / clustering / UMAP
# -----------------------------
pca_features <- VariableFeatures(obj)
if (length(pca_features) < 2) stop("Too few variable features for PCA.")

npcs_use <- min(as.integer(PCA_NPCS), length(pca_features))
npcs_use <- max(npcs_use, 2L)

obj <- RunPCA(
  obj,
  npcs = npcs_use,
  features = pca_features,
  verbose = FALSE
)

dims_nn   <- NEIGHBOR_DIMS[NEIGHBOR_DIMS <= npcs_use]
dims_umap <- UMAP_DIMS[UMAP_DIMS <= npcs_use]

if (length(dims_nn) < 2) stop("Too few PCA dimensions for FindNeighbors.")
if (length(dims_umap) < 2) stop("Too few PCA dimensions for RunUMAP.")

obj <- FindNeighbors(
  obj,
  dims = dims_nn,
  verbose = FALSE
)

obj <- FindClusters(
  obj,
  resolution = FINAL_RES,
  algorithm = 1,
  random.seed = 12345,
  verbose = FALSE
)

obj$seurat_clusters <- as.character(Idents(obj))

obj <- RunUMAP(
  obj,
  dims = dims_umap,
  n.neighbors = UMAP_N_NEIGHBORS,
  min.dist = UMAP_MIN_DIST,
  umap.method = "uwot",
  seed.use = 12345,
  verbose = FALSE
)

# -----------------------------
# Cluster labels
# -----------------------------
tcr_cluster_labels <- c(
  "C0" = "Mixed activated CD4/KLRB1+ state",
  "C1" = "Cycling T-cell state I",
  "C2" = "Cytokine-producing effector state",
  "C3" = "Cycling/effector T-cell state II",
  "C4" = "CD8/TRDC-associated cytotoxic state",
  "C5" = "CCR7/IL7R/HLA-II-associated state"
)

cluster_ids <- sort(unique(as.integer(as.character(obj$seurat_clusters))))
if (!identical(cluster_ids, 0:5)) {
  stop(
    "The frozen repeated-stimulation annotation requires clusters 0-5; observed: ",
    paste(cluster_ids, collapse = ", ")
  )
}
cluster_labels_short <- setNames(paste0("C", cluster_ids), as.character(cluster_ids))

expected_tcr_clusters <- paste0("C", cluster_ids)
missing_tcr_labels <- setdiff(expected_tcr_clusters, names(tcr_cluster_labels))
if (length(missing_tcr_labels) > 0) {
  stop("Missing TCR annotation(s): ", paste(missing_tcr_labels, collapse = ", "))
}

obj$cluster_short <- factor(
  unname(cluster_labels_short[as.character(obj$seurat_clusters)]),
  levels = paste0("C", cluster_ids)
)

obj$cluster_annot <- factor(
  unname(tcr_cluster_labels[as.character(obj$cluster_short)]),
  levels = unname(tcr_cluster_labels[paste0("C", cluster_ids)])
)

cluster_colors <- setNames(
  rep(cluster_palette, length.out = length(cluster_ids)),
  paste0("C", cluster_ids)
)

cluster_colors_annot <- setNames(
  unname(cluster_colors[paste0("C", cluster_ids)]),
  unname(tcr_cluster_labels[paste0("C", cluster_ids)])
)

# -----------------------------
# -----------------------------
# FIGURE 1: UMAP by clusters
# -----------------------------

# short labels only for text inside the UMAP
short_umap_labels <- c(
  "Mixed activated\nCD4/KLRB1+",
  "Cycling I",
  "Cytokine-producing\neffector",
  "Cycling/effector II",
  "CD8/TRDC-associated\ncytotoxic",
  "CCR7/IL7R/\nHLA-II-associated"
)

# Data-derived label positions avoid hard-coded coordinates and remain valid
# if the UMAP layout changes after a package update.
umap_coordinates <- as.data.frame(Embeddings(obj, "umap"))
colnames(umap_coordinates)[1:2] <- c("umap_1", "umap_2")
centers <- obj@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::bind_cols(umap_coordinates) %>%
  dplyr::group_by(cluster_short, cluster_annot) %>%
  dplyr::summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    short_label = short_umap_labels[match(as.character(cluster_short), expected_tcr_clusters)]
  )

p_umap_clusters <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "cluster_annot",
  cols = cluster_colors_annot,
  label = FALSE,
  pt.size = 0.5
) +
  geom_label(
    data = centers,
    aes(x = umap_1, y = umap_2, label = short_label),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 4.2,
    lineheight = 0.95,
    linewidth = 0,
    fill = scales::alpha("white", 0.75),
    label.padding = unit(0.10, "lines")
  ) +
  ggtitle("TCR UMAP by clusters") +
  theme_bw(base_size = 16) +
  theme(
    plot.title   = element_text(face = "bold", size = 20),
    axis.title   = element_text(face = "bold", size = 16),
    axis.text    = element_text(size = 13, color = "black"),
    legend.title = element_blank(),
    legend.text  = element_text(size = 13, face = "bold"),
    legend.key.height = unit(0.9, "cm")
  )

save_plot(p_umap_clusters, "Figure_5A_TCR_UMAP_clusters.png", w = 11, h = 8)

# -----------------------------
# FIGURE 2: cluster composition
# -----------------------------
md <- obj@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::mutate(
    sample = "TCR",
    cluster_short = factor(
      unname(cluster_labels_short[as.character(seurat_clusters)]),
      levels = paste0("C", cluster_ids)
    )
  )

df_bar_tcr <- md %>%
  dplyr::group_by(cluster_short) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(frac = n / sum(n))

df_bar_tcr$cluster_annot <- factor(
  unname(tcr_cluster_labels[as.character(df_bar_tcr$cluster_short)]),
  levels = unname(tcr_cluster_labels[paste0("C", cluster_ids)])
)

p_bar_tcr <- ggplot(df_bar_tcr, aes(x = cluster_annot, y = frac, fill = cluster_annot)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(values = cluster_colors_annot) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Repeated-stimulation cluster composition",
    x = "Cluster",
    y = "Fraction of QC-passing cells"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title   = element_text(face = "bold", size = 18),
    axis.title   = element_text(face = "bold", size = 16),
    axis.text    = element_text(size = 13, color = "black"),
    legend.position = "none"
  )

save_plot(p_bar_tcr, "Figure_5A_TCR_cluster_composition.png", w = 7, h = 5)

# -----------------------------
# FIGURE 3: QC violin plots
# -----------------------------
p_qc <- VlnPlot(
  obj,
  features = c("nFeature_RNA", "nCount_RNA"),
  ncol = 2,
  pt.size = 0
) &
  theme_bw(base_size = 14) &
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 14),
    axis.text  = element_text(size = 11, color = "black")
  )

save_plot(p_qc, "Figure_5A_TCR_QC_violin_optional.png", w = 10, h = 5)

# -----------------------------
# TABLE 1: Cell counts by cluster
# -----------------------------
cluster_cell_counts <- md %>%
  dplyr::group_by(cluster_short) %>%
  dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") %>%
  dplyr::arrange(cluster_short) %>%
  dplyr::mutate(fraction = n_cells / sum(n_cells))

write.csv(
  cluster_cell_counts,
  file.path(FIG_DIR, "Supplementary_Table_TCR_cluster_cell_counts.csv"),
  row.names = FALSE
)

# -----------------------------
# Exploratory cell-level marker ranking for cluster annotation. These tests do
# not constitute biological-replicate inference between experimental groups.
# -----------------------------
DefaultAssay(obj) <- "RNA"

markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

if (nrow(markers) > 0) {
  marker_counts <- markers %>%
    dplyr::mutate(cluster = as.character(.data$cluster)) %>%
    dplyr::count(.data$cluster, name = "n_markers")
  expected_marker_clusters <- as.character(0:5)
  if (!setequal(marker_counts$cluster, expected_marker_clusters) ||
      any(marker_counts$n_markers[match(expected_marker_clusters, marker_counts$cluster)] < 20L)) {
    stop("Every TCR cluster 0-5 must have at least 20 ranked markers for annotation review.")
  }
  fc_col <- if ("avg_log2FC" %in% colnames(markers)) {
    "avg_log2FC"
  } else if ("avg_logFC" %in% colnames(markers)) {
    "avg_logFC"
  } else {
    NULL
  }

  if (!is.null(fc_col)) {
    top_markers <- markers %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_max(order_by = .data[[fc_col]], n = 20, with_ties = FALSE) %>%
      dplyr::ungroup()

    wb <- openxlsx::createWorkbook()

    openxlsx::addWorksheet(wb, "All_clusters_top20")
    openxlsx::writeData(wb, "All_clusters_top20", top_markers)

    clusters <- sort(unique(top_markers$cluster))

    for (cl in clusters) {
      df_cl <- top_markers %>% dplyr::filter(cluster == cl)
      sheet_name <- paste0("Cluster_", cl)
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeData(wb, sheet_name, df_cl)
    }

    openxlsx::saveWorkbook(
      wb,
      file.path(FIG_DIR, "Supplementary_Table_TCR_top_markers_per_cluster.xlsx"),
      overwrite = TRUE
    )
  } else {
    stop("Could not find avg_log2FC or avg_logFC column in TCR marker table.")
  }
} else {
  stop("TCR marker table is empty; cluster annotations cannot be release-validated.")
}

# -----------------------------
# Save object
# -----------------------------
saveRDS(obj, file.path(FIG_DIR, "TCR_seurat_object.rds"))
writeLines(capture.output(sessionInfo()), file.path(FIG_DIR, "sessionInfo_R05.txt"))

# Both marker guards have now passed. Only validated outputs become visible in
# the stable publication directory.
promote_staged_outputs(FIG_DIR, FINAL_FIG_DIR)
unlink(FIG_DIR, recursive = TRUE, force = TRUE)

# -----------------------------
# Done
message("DONE.")
message("Source subset: ", SOURCE_SUBSET)
message("All outputs saved in: ", FINAL_FIG_DIR)
