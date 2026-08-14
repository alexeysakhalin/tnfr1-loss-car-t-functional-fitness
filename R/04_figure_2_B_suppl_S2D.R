
# ============================================================
# FIGURE 2B + SUPPLEMENTARY FIGURE S2D: KO vs WT RNA-seq analysis
#
# This script reproduces:
# - Figure 2B: Volcano plots for KO treated vs WT treated under:
#     • TNF
#     • IFNγ
#     • TNF + IFNγ
# - Supplementary Figure S2D: Venn diagram of downregulated genes
#   across the same three treatment conditions
#
# ------------------------------------------------------------
# DATA REQUIREMENTS (IMPORTANT)
#
# To run this script, obtain the provenance-documented RNA-seq inputs from the
# final public data release after its DOI is activated, or generate them from
# the authors' raw data using the declared pipeline and checksums.
#
# The RNA-seq dataset used here is located in:
#
#   data/rnaseq/
#
# Required file:
#   PyDESeq2_T6_vs_Hela_matched_treatments (1).xlsx
#
# Do NOT change file paths — the script assumes this structure.
#
# ------------------------------------------------------------
# OUTPUT
#
# The script generates:
# - Volcano plots for KO vs WT treated comparisons
# - Supplementary Venn diagram showing overlap of downregulated genes
#
# Output files are saved to:
#
#   figures/
#
# ------------------------------------------------------------
# NOTE
#
# Significant genes are defined as:
#   log2FC > 1 and adjusted p-value < 0.05   for upregulated genes
#   log2FC < -1 and adjusted p-value < 0.05  for downregulated genes
#
# Selected genes are highlighted on volcano plots.
# ============================================================

library(ggplot2)
library(dplyr)
library(readxl)
library(ggrepel)
library(grid)

dir.create("figures", recursive = TRUE, showWarnings = FALSE)

file_path <- file.path("data", "rnaseq", "PyDESeq2_T6_vs_Hela_matched_treatments (1).xlsx")

if (!file.exists(file_path)) {
  stop("Input file not found in data/rnaseq/")
}

make_volcano_G <- function(sheet_name, panel_title, out_name,
                           x_limits = c(-15, 15), y_max = 320) {

  genes <- read_excel(file_path, sheet = sheet_name) %>%
    filter(baseMean >= 30) %>%
    filter(!is.na(log2FoldChange), !is.na(padj)) %>%
    mutate(
      Gene_Symbol = toupper(trimws(Gene_Symbol)),
      padj_plot = ifelse(padj <= 0, 1e-300, padj)
    ) %>%
    mutate(
      status = case_when(
        log2FoldChange >  1 & padj < 0.05 ~ "Upregulated",
        log2FoldChange < -1 & padj < 0.05 ~ "Downregulated",
        TRUE                              ~ "Not significant"
      )
    ) %>%
    filter(is.finite(log2FoldChange), is.finite(-log10(padj_plot)))

  # same genes as Volcano A
  apoptosis_genes   <- c("CASP3","CASP7","CASP8","CASP9","BAX","BAK1","BCL2","FAS","APAF1")
  pyroptosis_genes  <- c("GSDMD","GSDME","CASP1","CASP4","CASP5","AIM2","NLRP3")
  necroptosis_genes <- c("RIPK1","RIPK3","MLKL")

  label_genes <- c(apoptosis_genes, pyroptosis_genes, necroptosis_genes, "ICAM1", "IRF1")

  label_data <- genes %>%
    filter(Gene_Symbol %in% label_genes) %>%
    mutate(
      y = -log10(padj_plot),
      side = ifelse(log2FoldChange < 0, "left", "right")
    ) %>%
    group_by(side) %>%
    arrange(y, .by_group = TRUE) %>%
    mutate(
      offset_rank = row_number(),
      nudge_x = ifelse(side == "left", -8.5, 8.5),
      nudge_y = (offset_rank - 1) * 8 + 6
    ) %>%
    ungroup()

  # panel-specific tuning
  if (sheet_name == "TNF") {
    pull_down <- c("CASP4","CASP7","CASP1","MLKL")

    label_data <- label_data %>%
      mutate(
        nudge_y = ifelse(Gene_Symbol == "ICAM1", pmax(0, nudge_y - 95), nudge_y),
        nudge_x = ifelse(Gene_Symbol == "ICAM1", nudge_x + 0.2, nudge_x),

        nudge_y = ifelse(Gene_Symbol == "CASP8", nudge_y + 8, nudge_y),
        nudge_x = ifelse(Gene_Symbol == "CASP8", nudge_x + 0.2, nudge_x),

        nudge_y = ifelse(Gene_Symbol == "IRF1", nudge_y + 18, nudge_y),
        nudge_x = ifelse(Gene_Symbol == "IRF1", nudge_x + 1.5, nudge_x),

        nudge_y = ifelse(Gene_Symbol %in% pull_down, pmax(0, nudge_y - 20), nudge_y)
      )
  }

  if (sheet_name == "IFN") {
    pull_down <- c("CASP4","CASP7","CASP1","MLKL")

    label_data <- label_data %>%
      mutate(
        nudge_y = ifelse(Gene_Symbol == "ICAM1", pmax(0, nudge_y - 55), nudge_y),
        nudge_x = ifelse(Gene_Symbol == "ICAM1", nudge_x + 1.5, nudge_x),

        nudge_y = ifelse(Gene_Symbol %in% pull_down, pmax(0, nudge_y - 40), nudge_y),

        nudge_y = ifelse(Gene_Symbol == "IRF1", pmax(0, nudge_y - 115), nudge_y),
        nudge_x = ifelse(Gene_Symbol == "IRF1", nudge_x + 1.2, nudge_x),

        nudge_x = ifelse(Gene_Symbol == "CASP4", nudge_x - 0.8, nudge_x),
        nudge_y = ifelse(Gene_Symbol == "CASP4", pmax(0, nudge_y - 10), nudge_y)
      )
  }

  if (sheet_name == "TI") {
    pull_down <- c("CASP4","CASP7","CASP1","MLKL")

    label_data <- label_data %>%
      mutate(
        nudge_y = ifelse(Gene_Symbol %in% pull_down, pmax(0, nudge_y - 40), nudge_y),
        nudge_y = ifelse(Gene_Symbol == "AIM2", nudge_y + 20, nudge_y),

        nudge_y = ifelse(Gene_Symbol == "ICAM1", pmax(0, nudge_y - 75), nudge_y),
        nudge_x = ifelse(Gene_Symbol == "ICAM1", nudge_x + 1.5, nudge_x),

        nudge_y = ifelse(Gene_Symbol == "IRF1", pmax(0, nudge_y - 125), nudge_y),
        nudge_x = ifelse(Gene_Symbol == "IRF1", nudge_x + 1.5, nudge_x),

        nudge_x = ifelse(Gene_Symbol == "CASP7", nudge_x + 0.8, nudge_x),
        nudge_y = ifelse(Gene_Symbol == "CASP7", nudge_y - 8, nudge_y)
      )
  }

  y_max_use <- max(
    y_max,
    max(label_data$y + label_data$nudge_y, na.rm = TRUE) + 5
  )

  p <- ggplot(genes, aes(x = log2FoldChange, y = -log10(padj_plot))) +
    geom_point(aes(color = status), size = 1.8, alpha = 0.85) +
    scale_color_manual(
      values = c(
        "Upregulated" = "#E64B35",
        "Downregulated" = "#3182BD",
        "Not significant" = "grey70"
      ),
      name = "Status"
    ) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    coord_cartesian(xlim = x_limits, ylim = c(0, y_max_use), clip = "on") +
    theme_minimal(base_size = 14) +
    labs(
      title = panel_title,
      x = "Log2 Fold Change",
      y = "-Log10 adjusted p-value"
    ) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "plain", size = 16),
      axis.title   = element_text(face = "plain", size = 14),
      axis.text    = element_text(size = 12),
      legend.title = element_text(face = "plain", size = 14),
      legend.text  = element_text(size = 12)
    ) +
    ggrepel::geom_label_repel(
      data = label_data,
      aes(label = Gene_Symbol),
      color = "black",
      size = 4,
      nudge_x = label_data$nudge_x,
      nudge_y = label_data$nudge_y,
      direction = "both",
      box.padding = 0.35,
      point.padding = 0.25,
      label.padding = grid::unit(0.15, "lines"),
      label.size = 0,
      fill = scales::alpha("white", 0.9),
      segment.color = "black",
      segment.size = 0.3,
      min.segment.length = 0,
      max.overlaps = Inf,
      force = 0.5,
      force_pull = 0.2,
      seed = 42,
      xlim = x_limits,
      ylim = c(0, y_max_use)
    )

  tiff(
    file.path("figures", paste0(out_name, ".tiff")),
    width = 8, height = 6, units = "in", res = 600
  )
  print(p)
  dev.off()

  png(
    file.path("figures", paste0(out_name, ".png")),
    width = 8, height = 6, units = "in", res = 600
  )
  print(p)
  dev.off()

  return(p)
}

# save 3 panels
p_ti <- make_volcano_G(
  sheet_name = "TI",
  panel_title = expression("TNF + IFN" * gamma),
  out_name = "Figure_2B_TNF_IFNg"
)

p_tnf <- make_volcano_G(
  sheet_name = "TNF",
  panel_title = "TNF",
  out_name = "Figure_2B_TNF"
)

p_ifn <- make_volcano_G(
  sheet_name = "IFN",
  panel_title = expression("IFN" * gamma),
  out_name = "Figure_2B_IFNg"
)














# =========================
# Supplementary S2D: Venn for DOWNREGULATED genes
# =========================

#install.packages(c("readxl", "dplyr", "VennDiagram", "writexl"))
library(readxl)
library(dplyr)
library(VennDiagram)
library(grid)
library(writexl)

file_path <- file.path("data", "rnaseq", "PyDESeq2_T6_vs_Hela_matched_treatments (1).xlsx")

# -------------------------
# 1. Load sheets
# -------------------------
deseq_tnf <- read_excel(file_path, sheet = "TNF")
deseq_ifn <- read_excel(file_path, sheet = "IFN")
deseq_ti  <- read_excel(file_path, sheet = "TI")

# -------------------------
# 2. Clean + select DOWN genes
# -------------------------
clean_tested_G <- function(df, gene_col = "Gene_Symbol") {
  df %>%
    mutate(
      Gene_Symbol = toupper(trimws(.data[[gene_col]]))
    ) %>%
    filter(baseMean >= 30) %>%
    filter(!is.na(Gene_Symbol), Gene_Symbol != "") %>%
    filter(!is.na(log2FoldChange), !is.na(padj)) %>%
    distinct(Gene_Symbol, .keep_all = TRUE)
}

clean_down_G <- function(df, gene_col = "Gene_Symbol") {
  clean_tested_G(df, gene_col) %>%
    filter(log2FoldChange < -1, padj < 0.05)
}

tnf_down <- clean_down_G(deseq_tnf)
ifn_down <- clean_down_G(deseq_ifn)
ti_down  <- clean_down_G(deseq_ti)

# -------------------------
# 3. Gene vectors
# -------------------------
down_tnf <- unique(tnf_down$Gene_Symbol)
down_ifn <- unique(ifn_down$Gene_Symbol)
down_ti  <- unique(ti_down$Gene_Symbol)

tested_universe_symbols <- Reduce(intersect, list(
  clean_tested_G(deseq_tnf)$Gene_Symbol,
  clean_tested_G(deseq_ifn)$Gene_Symbol,
  clean_tested_G(deseq_ti)$Gene_Symbol
))

length(down_tnf)
length(down_ifn)
length(down_ti)

length(intersect(down_tnf, down_ifn))
length(intersect(down_ifn, down_ti))
length(intersect(down_tnf, down_ti))
length(Reduce(intersect, list(down_tnf, down_ifn, down_ti)))

# -------------------------
# 4. Venn counts
# -------------------------
n12  <- length(intersect(down_tnf, down_ifn))
n23  <- length(intersect(down_ifn, down_ti))
n13  <- length(intersect(down_tnf, down_ti))
n123 <- length(Reduce(intersect, list(down_tnf, down_ifn, down_ti)))

# -------------------------
# 5. Draw Venn
# -------------------------
venn.plot <- draw.triple.venn(
  area1 = length(down_tnf),
  area2 = length(down_ifn),
  area3 = length(down_ti),
  n12 = n12,
  n23 = n23,
  n13 = n13,
  n123 = n123,
  category = c("TNF", "IFNγ", "TNF + IFNγ"),

  fill = c("#6FC7CF", "#1CC5FE", "#FBA27D"),
  alpha = 0.75,

  cex = 2,
  fontface = "bold",

  cat.cex = 2,
  cat.fontface = "bold",

  lwd = 2,
  scaled = FALSE
)

# -------------------------
# 6. Save Venn
# -------------------------

tiff(
  file.path("figures", "Supplementary_Figure_S2D_Venn_downregulated.tiff"),
  width = 6, height = 6, units = "in", res = 600
)
grid.draw(venn.plot)
dev.off()

png(
  file.path("figures", "Supplementary_Figure_S2D_Venn_downregulated.png"),
  width = 6, height = 6, units = "in", res = 600
)
grid.draw(venn.plot)
dev.off()












# =========================
# Suppl S2D making Excel file for common DOWN genes
# =========================

common_down <- Reduce(intersect, list(down_tnf, down_ifn, down_ti))
length(common_down)

common_down_genes_only <- data.frame(Gene_Symbol = sort(common_down))

tnf_tbl <- tnf_down %>%
  dplyr::select(Gene_Symbol, log2FoldChange, padj, baseMean) %>%
  dplyr::rename(
    log2FC_TNF = log2FoldChange,
    padj_TNF = padj,
    baseMean_TNF = baseMean
  )

ifn_tbl <- ifn_down %>%
  dplyr::select(Gene_Symbol, log2FoldChange, padj, baseMean) %>%
  dplyr::rename(
    log2FC_IFNg = log2FoldChange,
    padj_IFNg = padj,
    baseMean_IFNg = baseMean
  )

ti_tbl <- ti_down %>%
  dplyr::select(Gene_Symbol, log2FoldChange, padj, baseMean) %>%
  dplyr::rename(
    log2FC_TNF_IFNg = log2FoldChange,
    padj_TNF_IFNg = padj,
    baseMean_TNF_IFNg = baseMean
  )

common_down_with_stats <- data.frame(Gene_Symbol = common_down) %>%
  left_join(tnf_tbl, by = "Gene_Symbol") %>%
  left_join(ifn_tbl, by = "Gene_Symbol") %>%
  left_join(ti_tbl, by = "Gene_Symbol") %>%
  arrange(log2FC_TNF_IFNg)

# -------------------------
# Save Excel (downregulated genes)
# -------------------------
save_path <- file.path("figures", "Supplementary_Figure_S2D_common_downregulated_genes.xlsx")

write_xlsx(
  list(
    "common_down_genes_only" = common_down_genes_only,
    "common_down_with_log2FC" = common_down_with_stats
  ),
  path = save_path
)

cat("Excel saved:\n")
cat(save_path, "\n")












# =========================
# Suppl S2D
# =========================
# =========================
# Enrichment: GO BP / CC / MF / KEGG
# for common downregulated genes (Optional)
# =========================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(openxlsx)
library(ggplot2)
library(dplyr)

genes_down <- common_down

# -------------------------
# 1. SYMBOL -> ENTREZ
# -------------------------
gene_df <- bitr(
  genes_down,
  fromType = "SYMBOL",
  toType = c("ENTREZID"),
  OrgDb = org.Hs.eg.db
)

entrez_genes <- unique(gene_df$ENTREZID)

universe_df <- bitr(
  tested_universe_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
entrez_universe <- unique(universe_df$ENTREZID)

# -------------------------
# 2. Enrichment
# -------------------------
ego_bp <- enrichGO(
  gene          = entrez_genes,
  universe      = entrez_universe,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.1,
  qvalueCutoff  = 1,
  readable      = TRUE
)

ego_cc <- enrichGO(
  gene          = entrez_genes,
  universe      = entrez_universe,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.1,
  qvalueCutoff  = 1,
  readable      = TRUE
)

ego_mf <- enrichGO(
  gene          = entrez_genes,
  universe      = entrez_universe,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.1,
  qvalueCutoff  = 1,
  readable      = TRUE
)

ego_kegg <- enrichKEGG(
  gene          = entrez_genes,
  universe      = entrez_universe,
  organism      = "hsa",
  pvalueCutoff  = 0.1,
  pAdjustMethod = "BH",
  qvalueCutoff  = 1
)

# -------------------------
# 3. Convert to data frames
# -------------------------
bp_df   <- as.data.frame(ego_bp)
cc_df   <- as.data.frame(ego_cc)
mf_df   <- as.data.frame(ego_mf)
kegg_df <- as.data.frame(ego_kegg)

cat("GO BP rows:", nrow(bp_df), "\n")
cat("GO CC rows:", nrow(cc_df), "\n")
cat("GO MF rows:", nrow(mf_df), "\n")
cat("KEGG rows:", nrow(kegg_df), "\n")

# -------------------------
# 4. Save all non-empty results to Excel
# -------------------------
out_list <- list(
  mapped_genes = gene_df
)

if (nrow(bp_df) > 0)   out_list[["GO_BP"]] <- bp_df
if (nrow(cc_df) > 0)   out_list[["GO_CC"]] <- cc_df
if (nrow(mf_df) > 0)   out_list[["GO_MF"]] <- mf_df
if (nrow(kegg_df) > 0) out_list[["KEGG"]]  <- kegg_df

# =========================
# 4. Save Excel
# =========================
excel_path <- file.path("figures", "Supplementary_Figure_S2D_common_down_enrichment_all.xlsx")

write.xlsx(
  out_list,
  file = excel_path,
  overwrite = TRUE
)

cat("Excel saved:\n")
cat(excel_path, "\n")


# =========================
# 5. Function to save dotplot + gene ratio plot
# =========================
save_enrichment_plots <- function(enrich_obj, enrich_df, prefix, title_text) {

  if (nrow(enrich_df) == 0) {
    cat(prefix, ": no enriched terms, skipping plots\n")
    return(NULL)
  }

  n_show <- min(12, nrow(enrich_df))

  # ---- dotplot ----
  p_dot <- dotplot(
    enrich_obj,
    showCategory = n_show,
    font.size = 12,
    title = title_text
  ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.y = element_text(size = 11),
      axis.text.x = element_text(size = 11),
      axis.title = element_text(size = 12)
    )

  tiff(
    file.path("figures", paste0(prefix, "_dotplot.tiff")),
    width = 9, height = 7, units = "in", res = 600
  )
  print(p_dot)
  dev.off()

  png(
    file.path("figures", paste0(prefix, "_dotplot.png")),
    width = 9, height = 7, units = "in", res = 600
  )
  print(p_dot)
  dev.off()

  # ---- gene ratio plot ----
  plot_df <- enrich_df %>%
    dplyr::slice_head(n = n_show) %>%
    dplyr::mutate(
      Description = factor(Description, levels = rev(Description)),
      GeneRatio_num = sapply(GeneRatio, function(x) {
        parts <- strsplit(x, "/")[[1]]
        as.numeric(parts[1]) / as.numeric(parts[2])
      })
    )

  p_ratio <- ggplot(plot_df, aes(x = GeneRatio_num, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    theme_bw(base_size = 12) +
    labs(
      title = paste0(title_text, " - Gene Ratio"),
      x = "Gene Ratio",
      y = NULL,
      color = "Adjusted p-value",
      size = "Gene count"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.y = element_text(size = 11),
      axis.text.x = element_text(size = 11),
      axis.title = element_text(size = 12)
    )

  tiff(
    file.path("figures", paste0(prefix, "_gene_ratio.tiff")),
    width = 9, height = 7, units = "in", res = 600
  )
  print(p_ratio)
  dev.off()

  png(
    file.path("figures", paste0(prefix, "_gene_ratio.png")),
    width = 9, height = 7, units = "in", res = 600
  )
  print(p_ratio)
  dev.off()
}

# =========================
# 6. Save plots for each category
# =========================
save_enrichment_plots(
  ego_bp, bp_df,
  prefix = "Supplementary_Figure_S2D_common_down_GO_BP",
  title_text = "GO Biological Process enrichment"
)

save_enrichment_plots(
  ego_cc, cc_df,
  prefix = "Supplementary_Figure_S2D_common_down_GO_CC",
  title_text = "GO Cellular Component enrichment"
)

save_enrichment_plots(
  ego_mf, mf_df,
  prefix = "Supplementary_Figure_S2D_common_down_GO_MF",
  title_text = "GO Molecular Function enrichment"
)

save_enrichment_plots(
  ego_kegg, kegg_df,
  prefix = "Supplementary_Figure_S2D_common_down_KEGG",
  title_text = "KEGG pathway enrichment"
)
