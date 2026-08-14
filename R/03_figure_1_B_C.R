
# ============================================================
# FIGURE 1B–C: RNA-seq differential expression analysis
#
# This script reproduces:
# - Figure 1B: Volcano plots (TNF, IFNγ, TNF+IFNγ)
# - Figure 1C: Overlap of upregulated genes (Venn diagram)
#
# ------------------------------------------------------------
# DATA REQUIREMENTS (IMPORTANT)
#
# To run this script, obtain the provenance-documented RNA-seq inputs from the
# final public data release after its DOI is activated, or generate them from
# the authors' raw data using the declared pipeline and checksums.
#
# The RNA-seq differential expression data used here are located in:
#
#   data/rnaseq/
#
# Required files:
#   Differential_Expression_TNFa_vs_control_final_filtered.xlsx; Differential_Expression_IFN_vs_control.xlsx;
# Differential_Expression_TI_vs_control_final_filtered.xlsx
#
# Do NOT change file paths — the script assumes this structure.
#
# ------------------------------------------------------------
# OUTPUT
#
# The script generates:
# - Volcano plots comparing treatment vs untreated conditions:
#     • TNF
#     • IFNγ
#     • TNF + IFNγ
#
# - Overlap of upregulated genes across conditions
#   (used to identify shared genes such as ICAM1 and IRF1)
#
# Corresponding to:
#   • Figure 1B (volcano plots)
#   • Figure 1C (gene overlap / Venn diagram)
#
# ------------------------------------------------------------
# NOTE
#
# Upregulated genes are defined as:
#   log2FC > 1 and adjusted p-value < 0.05
#
# Gene labels in volcano plots highlight key apoptosis
# and inflammatory regulators.
# ============================================================

library(ggplot2)
library(dplyr)
library(readxl)
library(ggrepel)

dir.create("figures", recursive = TRUE, showWarnings = FALSE)

genes <- readxl::read_excel(
  file.path("data", "rnaseq", "Differential_Expression_TNFa_vs_control_final_filtered.xlsx")
) %>%
  filter(baseMean >= 30) %>%
  filter(!is.na(log2FoldChange), !is.na(padj)) %>%
  mutate(Gene_ID = toupper(trimws(Gene_ID)))

genes <- genes %>%
  mutate(status = case_when(
    log2FoldChange >  1 & padj < 0.05 ~ "Upregulated",
    log2FoldChange < -1 & padj < 0.05 ~ "Downregulated",
    TRUE                               ~ "Not significant"
  ))

genes %>% filter(Gene_ID == "ICAM1")

apoptosis_genes  <- c("CASP3","CASP7","CASP8","CASP9","BAX","BAK1","BCL2","FAS","APAF1")
pyroptosis_genes <- c("GSDMD","GSDME","CASP1","CASP4","CASP5","AIM2","NLRP3")
necroptosis_genes<- c("RIPK1","RIPK3","MLKL")

label_genes <- c(apoptosis_genes, pyroptosis_genes, necroptosis_genes, "ICAM1", "IRF1")
genes %>% filter(Gene_ID == "IRF1")

dx <- 9    #
dy <- 9    #

label_data <- genes %>%
  filter(Gene_ID %in% label_genes) %>%
  mutate(
    y    = -log10(pmax(padj, .Machine$double.xmin)),
    side = ifelse(log2FoldChange < 0, "left", "right")
  ) %>%
  group_by(side) %>%
  arrange(y, .by_group = TRUE) %>%
  mutate(
    offset_rank = row_number(),
    nudge_x = ifelse(side == "left", -dx, dx),   #
    nudge_y = (offset_rank - 1) * dy + 6         #
  ) %>%
  ungroup()

label_data <- label_data %>%
  mutate(
    nudge_y = ifelse(Gene_ID == "ICAM1", nudge_y - 95, nudge_y),
    nudge_x = ifelse(Gene_ID == "ICAM1", nudge_x + 0.2, nudge_x)
  )

label_data <- label_data %>%
  mutate(
    nudge_y = ifelse(Gene_ID == "IRF1", nudge_y + 12, nudge_y),
    nudge_x = ifelse(Gene_ID == "IRF1", nudge_x + 1.0, nudge_x)
  )

y_max <- 300
x_limits <- c(-15, 15)

label_data <- label_data %>%
  dplyr::mutate(
    nudge_y = ifelse(Gene_ID == "CASP8", nudge_y + 8, nudge_y),
    nudge_x = ifelse(Gene_ID == "CASP8", nudge_x + 0.2, nudge_x) #
  )

y_max <- max(
  stats::quantile(-log10(genes$padj), 0.99, na.rm = TRUE),
  max(label_data$y + label_data$nudge_y, na.rm = TRUE)
) + 5

volcano <- ggplot(genes, aes(x = log2FoldChange,
                            y = -log10(pmax(padj, .Machine$double.xmin)))) +
  geom_point(aes(color = status), size = 1.8, alpha = 0.85) +
  scale_color_manual(values = c("Upregulated" = "#E64B35",
                                "Downregulated" = "#3182BD",
                                "Not significant" = "grey70"),
                     name = "Status") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  coord_cartesian(xlim = x_limits, ylim = c(0, y_max), clip = "on") +
  theme_minimal(base_size = 14) +
  labs(title = "TNF",
       x = "Log2 Fold Change",
       y = "-Log10 adjusted p-value") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain")
  ) +
  ggrepel::geom_label_repel(
    data = label_data,
    aes(label = Gene_ID),
    color = "black", size = 4,
    nudge_x = label_data$nudge_x,
    nudge_y = label_data$nudge_y,
    direction = "both",
    box.padding = 0.35, point.padding = 0.25,   #
    label.padding = grid::unit(0.15, "lines"),
    label.size = 0,                              #
    fill = scales::alpha("white", 0.9),          #
    segment.color = "black", segment.size = 0.3,
    max.overlaps = Inf, force = 0.5, force_pull = 0.2,
    seed = 42,
    xlim = x_limits, ylim = c(0, y_max)
  )

# ---------- FIGURE 1B: TNF ----------
tiff("figures/Figure_1B_TNF_volcano.tiff",
     width = 8, height = 6, units = "in", res = 600)
print(volcano)
dev.off()

png("figures/Figure_1B_TNF_volcano.png",
    width = 8, height = 6, units = "in", res = 600)
print(volcano)
dev.off()






#IFN

library(ggplot2)
library(dplyr)
library(readxl)
library(ggrepel)

ifn_file <- file.path(
  "data", "rnaseq",
  "Differential_Expression_IFN_vs_control.xlsx"
)

if (!file.exists(ifn_file)) {
  stop("IFN RNA-seq file not found in data/rnaseq/")
}

genes <- readxl::read_excel(ifn_file) %>%
  dplyr::filter(baseMean >= 30) %>%
  dplyr::filter(!is.na(log2FoldChange), !is.na(padj)) %>%
  dplyr::mutate(Gene_ID = toupper(trimws(Gene_ID)))

genes <- genes %>%
  mutate(status = case_when(
    log2FoldChange >  1 & padj < 0.05 ~ "Upregulated",
    log2FoldChange < -1 & padj < 0.05 ~ "Downregulated",
    TRUE                               ~ "Not significant"
  ))

genes %>% filter(Gene_ID == "ICAM1")

apoptosis_genes  <- c("CASP3","CASP7","CASP8","CASP9","BAX","BAK1","BCL2","FAS","APAF1")
pyroptosis_genes <- c("GSDMD","GSDME","CASP1","CASP4","CASP5","AIM2","NLRP3")
necroptosis_genes<- c("RIPK1","RIPK3","MLKL")

label_genes <- c(apoptosis_genes, pyroptosis_genes, necroptosis_genes, "ICAM1", "IRF1")
genes %>% filter(Gene_ID == "IRF1")

dx <- 9   #
dy <- 9   #

label_data <- genes %>%
  filter(Gene_ID %in% label_genes) %>%
  mutate(
    y    = -log10(pmax(padj, .Machine$double.xmin)),
    side = ifelse(log2FoldChange < 0, "left", "right")
  ) %>%
  group_by(side) %>%
  arrange(y, .by_group = TRUE) %>%
  mutate(
    offset_rank = row_number(),
    nudge_x = ifelse(side == "left", -dx, dx),
    nudge_y = (offset_rank - 1) * dy + 6
  ) %>%
  ungroup()

label_data <- label_data %>%
  mutate(
    nudge_y = ifelse(Gene_ID == "ICAM1", pmax(0, nudge_y - 55), nudge_y),
    nudge_x = ifelse(Gene_ID == "ICAM1", nudge_x + 1.5, nudge_x)
  )

label_data <- label_data %>%
  mutate(
    nudge_y = ifelse(Gene_ID == "IRF1", pmax(0, nudge_y - 115), nudge_y),
    nudge_x = ifelse(Gene_ID == "IRF1", nudge_x + 1.2, nudge_x)
  )

pull_down <- c("CASP4","CASP7","CASP1","MLKL")
label_data <- label_data %>%
  mutate(
    nudge_y = ifelse(Gene_ID %in% pull_down, pmax(0, nudge_y - 40), nudge_y) #
  )

y_max <- 300
x_limits <- c(-15, 15)

volcano <- ggplot(genes, aes(x = log2FoldChange,
                            y = -log10(pmax(padj, .Machine$double.xmin)))) +
  geom_point(aes(color = status), size = 1.8, alpha = 0.75) +   #
  scale_color_manual(values = c("Upregulated" = "#E64B35",
                                "Downregulated" = "#3182BD",
                                "Not significant" = "grey70"),
                     name = "Status") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  coord_cartesian(xlim = x_limits, ylim = c(0, y_max), clip = "on") +
  theme_minimal(base_size = 14) +
  labs(title = expression("IFN" * gamma),
       x = "Log2 Fold Change",
       y = "-Log10 adjusted p-value") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain")
  ) +
  ggrepel::geom_label_repel(
    data = label_data,
    aes(label = Gene_ID),
    color = "black", size = 4,
    nudge_x = label_data$nudge_x,
    nudge_y = label_data$nudge_y,
    direction = "both",
    box.padding = 0.35,
    point.padding = 0.6,
    label.padding = grid::unit(0.15, "lines"),
    label.size = 0,
    fill = scales::alpha("white", 0.9),
    segment.color = "black",
    segment.size = 0.3,
    min.segment.length = 0,
    max.overlaps = Inf,
    force = 0.6,
    force_pull = 0.2,
    seed = 42,
    xlim = x_limits, ylim = c(0, y_max)
  )

# ---------- FIGURE 1B: IFNγ ----------
tiff("figures/Figure_1B_IFNg_volcano.tiff",
     width = 8, height = 6, units = "in", res = 600)
print(volcano)
dev.off()

png("figures/Figure_1B_IFNg_volcano.png",
    width = 8, height = 6, units = "in", res = 600)
print(volcano)
dev.off()









#TNF + IFN

library(ggplot2)
library(dplyr)
library(readxl)
library(ggrepel)

ti_file <- file.path(
  "data", "rnaseq",
  "Differential_Expression_TI_vs_control_final_filtered.xlsx"
)

if (!file.exists(ti_file)) {
  stop("TNF+IFNγ RNA-seq file not found in data/rnaseq/")
}

genes <- readxl::read_excel(ti_file) %>%
  dplyr::filter(baseMean >= 30) %>%
  dplyr::filter(!is.na(log2FoldChange), !is.na(padj)) %>%
  dplyr::mutate(Gene_ID = toupper(trimws(Gene_ID)))

genes <- genes %>%
  mutate(status = case_when(
    log2FoldChange >  1 & padj < 0.05 ~ "Upregulated",
    log2FoldChange < -1 & padj < 0.05 ~ "Downregulated",
    TRUE                               ~ "Not significant"
  ))

genes %>% filter(Gene_ID == "ICAM1")

apoptosis_genes  <- c("CASP3","CASP7","CASP8","CASP9","BAX","BAK1","BCL2","FAS","APAF1")
pyroptosis_genes <- c("GSDMD","GSDME","CASP1","CASP4","CASP5","AIM2","NLRP3")
necroptosis_genes<- c("RIPK1","RIPK3","MLKL")

label_genes <- c(apoptosis_genes, pyroptosis_genes, necroptosis_genes, "ICAM1", "IRF1")
genes %>% filter(Gene_ID == "IRF1")

dx <- 9   #
dy <- 9   #

label_data <- genes %>%
  filter(Gene_ID %in% label_genes) %>%
  mutate(
    y    = -log10(pmax(padj, .Machine$double.xmin)),
    side = ifelse(log2FoldChange < 0, "left", "right")
  ) %>%
  group_by(side) %>%
  arrange(y, .by_group = TRUE) %>%
  mutate(
    offset_rank = row_number(),
    nudge_x = ifelse(side == "left", -dx, dx),
    nudge_y = (offset_rank - 1) * dy + 6
  ) %>%
  ungroup()

label_data <- label_data %>%
  mutate(
    nudge_y = ifelse(Gene_ID %in% pull_down, pmax(0, nudge_y - 40), nudge_y),
    nudge_y = ifelse(Gene_ID == "AIM2", nudge_y + 20, nudge_y),
    nudge_y = ifelse(Gene_ID == "ICAM1", pmax(0, nudge_y - 75), nudge_y),
    nudge_x = ifelse(Gene_ID == "ICAM1", nudge_x + 1.5, nudge_x)
  )

label_data <- label_data %>%
  mutate(
    nudge_y = ifelse(Gene_ID == "IRF1", pmax(0, nudge_y - 125), nudge_y),
    nudge_x = ifelse(Gene_ID == "IRF1", nudge_x + 1.5, nudge_x)
  )

pull_down <- c("CASP4","CASP7","CASP1","MLKL")
label_data <- label_data %>%
  mutate(
    nudge_y = ifelse(Gene_ID %in% pull_down, pmax(0, nudge_y - 40), nudge_y) #
  )

y_max <- 300
x_limits <- c(-15, 15)

label_data <- label_data %>%
  mutate(
    nudge_y = ifelse(Gene_ID %in% pull_down, pmax(0, nudge_y - 40), nudge_y),
    nudge_y = ifelse(Gene_ID == "AIM2", nudge_y + 20, nudge_y)  #
  )

label_data <- label_data %>%
  mutate(
    nudge_x = ifelse(Gene_ID == "CASP7", nudge_x + 0.8, nudge_x),
    nudge_y = ifelse(Gene_ID == "CASP7", nudge_y - 8,  nudge_y)
  )

volcano <- ggplot(genes, aes(x = log2FoldChange,
                            y = -log10(pmax(padj, .Machine$double.xmin)))) +
  geom_point(aes(color = status), size = 1.8, alpha = 0.75) +   #
  scale_color_manual(values = c("Upregulated" = "#E64B35",
                                "Downregulated" = "#3182BD",
                                "Not significant" = "grey70"),
                     name = "Status") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  coord_cartesian(xlim = x_limits, ylim = c(0, y_max), clip = "on") +
  theme_minimal(base_size = 14) +
  labs(title = expression("TNF + IFN" * gamma),
       x = "Log2 Fold Change",
       y = "-Log10 adjusted p-value") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain")
  ) +
  ggrepel::geom_label_repel(
    data = label_data,
    aes(label = Gene_ID),
    color = "black", size = 4,
    nudge_x = label_data$nudge_x,
    nudge_y = label_data$nudge_y,
    direction = "both",
    box.padding = 0.35, point.padding = 0.6,
    label.padding = grid::unit(0.15, "lines"),
    label.size = 0, fill = scales::alpha("white", 0.9),
    segment.color = "black", segment.size = 0.3,
    min.segment.length = 0,
    max.overlaps = Inf, force = 0.6, force_pull = 0.2,
    seed = 42,
    xlim = x_limits, ylim = c(0, y_max)
  )

# ---------- FIGURE 1B: TNF + IFNγ ----------
tiff("figures/Figure_1B_TNF_IFNg_volcano.tiff",
     width = 8, height = 6, units = "in", res = 600)
print(volcano)
dev.off()

png("figures/Figure_1B_TNF_IFNg_volcano.png",
    width = 8, height = 6, units = "in", res = 600)
print(volcano)
dev.off()









# Figure_1_C_Venn_diagram

library(readxl)
library(dplyr)
library(writexl)
library(VennDiagram)
library(grid)

deseq_tnf <- read_excel(
  file.path("data", "rnaseq", "Differential_Expression_TNFa_vs_control_final_filtered.xlsx")
)

deseq_ifng <- read_excel(
  file.path("data", "rnaseq", "Differential_Expression_IFN_vs_control.xlsx")
)

deseq_ti <- read_excel(
  file.path("data", "rnaseq", "Differential_Expression_TI_vs_control_final_filtered.xlsx")
)

colnames(deseq_tnf)
colnames(deseq_ifng)
colnames(deseq_ti)

deseq_tnf  <- deseq_tnf  %>% filter(baseMean >= 30)
deseq_ifng <- deseq_ifng %>% filter(baseMean >= 30)
deseq_ti   <- deseq_ti   %>% filter(baseMean >= 30)

head(deseq_tnf)
head(deseq_ifng)
head(deseq_ti)

sum(deseq_tnf$log2FoldChange > 1 & deseq_tnf$padj < 0.05, na.rm = TRUE)
sum(deseq_ifng$log2FoldChange > 1 & deseq_ifng$padj < 0.05, na.rm = TRUE)
sum(deseq_ti$log2FoldChange > 1 & deseq_ti$padj < 0.05, na.rm = TRUE)

up_tnf <- deseq_tnf$Gene_ID[deseq_tnf$log2FoldChange > 1 & deseq_tnf$padj < 0.05]
up_ifng <- deseq_ifng$Gene_ID[deseq_ifng$log2FoldChange > 1 & deseq_ifng$padj < 0.05]
up_ti <- deseq_ti$Gene_ID[deseq_ti$log2FoldChange > 1 & deseq_ti$padj < 0.05]

up_tnf  <- unique(toupper(trimws(up_tnf)))
up_ifng <- unique(toupper(trimws(up_ifng)))
up_ti   <- unique(toupper(trimws(up_ti)))

up_tnf  <- up_tnf[!is.na(up_tnf) & up_tnf != ""]
up_ifng <- up_ifng[!is.na(up_ifng) & up_ifng != ""]
up_ti   <- up_ti[!is.na(up_ti) & up_ti != ""]

length(up_tnf)
length(up_ifng)
length(up_ti)

length(intersect(up_tnf, up_ifng))
length(intersect(up_ifng, up_ti))
length(intersect(up_tnf, up_ti))
length(Reduce(intersect, list(up_tnf, up_ifng, up_ti)))

n12 <- length(intersect(up_tnf, up_ifng))
n23 <- length(intersect(up_ifng, up_ti))
n13 <- length(intersect(up_tnf, up_ti))
n123 <- length(Reduce(intersect, list(up_tnf, up_ifng, up_ti)))

venn.plot <- draw.triple.venn(
  area1 = length(up_tnf),
  area2 = length(up_ifng),
  area3 = length(up_ti),
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
  cat.fontface = "bold"
)

tiff(
  filename = file.path("figures", "Figure_1_C.tiff"),
  width = 6, height = 6, units = "in", res = 600
)
grid.draw(venn.plot)
dev.off()

png(
  filename = file.path("figures", "Figure_1_C.png"),
  width = 6, height = 6, units = "in", res = 600
)
grid.draw(venn.plot)
dev.off()










# making Excel file for 187 shared upregulated genes (Optional)

# install.packages(c("readxl", "dplyr", "writexl"))
library(readxl)
library(dplyr)
library(writexl)

deseq_tnf <- read_excel(
  file.path("data", "rnaseq", "Differential_Expression_TNFa_vs_control_final_filtered.xlsx")
)

deseq_ifng <- read_excel(
  file.path("data", "rnaseq", "Differential_Expression_IFN_vs_control.xlsx")
)

deseq_ti <- read_excel(
  file.path("data", "rnaseq", "Differential_Expression_TI_vs_control_final_filtered.xlsx")
)

clean_tested <- function(df, gene_col = "Gene_ID") {
  df %>%
    mutate(
      Gene_ID = toupper(trimws(.data[[gene_col]]))
    ) %>%
    filter(baseMean >= 30) %>%
    filter(!is.na(Gene_ID), Gene_ID != "") %>%
    filter(!is.na(log2FoldChange), !is.na(padj)) %>%
    distinct(Gene_ID, .keep_all = TRUE)
}

clean_up <- function(df, gene_col = "Gene_ID") {
  clean_tested(df, gene_col) %>%
    filter(log2FoldChange > 1, padj < 0.05)
}

tnf_up <- clean_up(deseq_tnf)
ifn_up <- clean_up(deseq_ifng)
ti_up  <- clean_up(deseq_ti)

length(unique(tnf_up$Gene_ID))
length(unique(ifn_up$Gene_ID))
length(unique(ti_up$Gene_ID))

common_187 <- Reduce(intersect, list(tnf_up$Gene_ID, ifn_up$Gene_ID, ti_up$Gene_ID))
length(common_187)   #

# Enrichment background: genes that were eligible/tested in all three DE
# contrasts, not the whole genome.
tested_universe_symbols <- Reduce(intersect, list(
  clean_tested(deseq_tnf)$Gene_ID,
  clean_tested(deseq_ifng)$Gene_ID,
  clean_tested(deseq_ti)$Gene_ID
))

common_genes_only <- data.frame(Gene_ID = sort(common_187))

tnf_tbl <- tnf_up %>%
  select(Gene_ID, log2FoldChange, padj, baseMean) %>%
  rename(
    log2FC_TNF = log2FoldChange,
    padj_TNF = padj,
    baseMean_TNF = baseMean
  )

ifn_tbl <- ifn_up %>%
  select(Gene_ID, log2FoldChange, padj, baseMean) %>%
  rename(
    log2FC_IFNg = log2FoldChange,
    padj_IFNg = padj,
    baseMean_IFNg = baseMean
  )

ti_tbl <- ti_up %>%
  select(Gene_ID, log2FoldChange, padj, baseMean) %>%
  rename(
    log2FC_TNF_IFNg = log2FoldChange,
    padj_TNF_IFNg = padj,
    baseMean_TNF_IFNg = baseMean
  )

common_with_stats <- data.frame(Gene_ID = common_187) %>%
  left_join(tnf_tbl, by = "Gene_ID") %>%
  left_join(ifn_tbl, by = "Gene_ID") %>%
  left_join(ti_tbl, by = "Gene_ID") %>%
  arrange(desc(log2FC_TNF_IFNg))

# =========================

save_path <- file.path("figures", "Figure_1_C_Common_187_genes_upregulated.xlsx")

write_xlsx(
  list(
    "common_187_genes_only" = common_genes_only,
    "common_187_with_log2FC" = common_with_stats
  ),
  path = save_path
)








# =========================
# Enrichment for
# common_187 upregulated genes (Optional)
# =========================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(openxlsx)
library(ggplot2)
library(dplyr)

genes_up <- common_187

# -------------------------
# 1. SYMBOL -> ENTREZ
# -------------------------
gene_df_up <- bitr(
  genes_up,
  fromType = "SYMBOL",
  toType = c("ENTREZID"),
  OrgDb = org.Hs.eg.db
)

entrez_up <- unique(gene_df_up$ENTREZID)

universe_df_up <- bitr(
  tested_universe_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
entrez_universe_up <- unique(universe_df_up$ENTREZID)

cat("Input genes:", length(genes_up), "\n")
cat("Mapped ENTREZ:", length(entrez_up), "\n")

# -------------------------
# 2. Enrichment
# -------------------------
ego_bp_up <- enrichGO(
  gene          = entrez_up,
  universe      = entrez_universe_up,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.1,
  qvalueCutoff  = 1,
  readable      = TRUE
)

ego_cc_up <- enrichGO(
  gene          = entrez_up,
  universe      = entrez_universe_up,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.1,
  qvalueCutoff  = 1,
  readable      = TRUE
)

ego_mf_up <- enrichGO(
  gene          = entrez_up,
  universe      = entrez_universe_up,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.1,
  qvalueCutoff  = 1,
  readable      = TRUE
)

ego_kegg_up <- enrichKEGG(
  gene          = entrez_up,
  universe      = entrez_universe_up,
  organism      = "hsa",
  pvalueCutoff  = 0.1,
  pAdjustMethod = "BH",
  qvalueCutoff  = 1
)

# -------------------------
# 3. Convert to data frames
# -------------------------
bp_up_df   <- as.data.frame(ego_bp_up)
cc_up_df   <- as.data.frame(ego_cc_up)
mf_up_df   <- as.data.frame(ego_mf_up)
kegg_up_df <- as.data.frame(ego_kegg_up)

cat("GO BP rows:", nrow(bp_up_df), "\n")
cat("GO CC rows:", nrow(cc_up_df), "\n")
cat("GO MF rows:", nrow(mf_up_df), "\n")
cat("KEGG rows:", nrow(kegg_up_df), "\n")

# -------------------------
# 4. Save Excel
# -------------------------
out_list_up <- list(
  mapped_genes = gene_df_up
)

if (nrow(bp_up_df) > 0)   out_list_up[["GO_BP"]] <- bp_up_df
if (nrow(cc_up_df) > 0)   out_list_up[["GO_CC"]] <- cc_up_df
if (nrow(mf_up_df) > 0)   out_list_up[["GO_MF"]] <- mf_up_df
if (nrow(kegg_up_df) > 0) out_list_up[["KEGG"]]  <- kegg_up_df

# =========================
excel_path <- file.path("figures", "Figure_1_C_Common_187_upregulated_enrichment_all.xlsx")

write.xlsx(
  out_list_up,
  file = excel_path,
  overwrite = TRUE
)

# =========================
# 5. Function to save dotplot + gene ratio plot
save_enrichment_plots <- function(enrich_obj, enrich_df, prefix, title_text, n_show = 12) {

  save_dir <- "figures"
  dir.create(save_dir, showWarnings = FALSE)

  if (nrow(enrich_df) == 0) {
    cat(prefix, ": no enriched terms, skipping plots\n")
    return(NULL)
  }

  n_show <- min(n_show, nrow(enrich_df))

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
    file.path(save_dir, paste0(prefix, "_dotplot.tiff")),
    width = 9, height = 7, units = "in", res = 600
  )
  print(p_dot)
  dev.off()

  png(
    file.path(save_dir, paste0(prefix, "_dotplot.png")),
    width = 9, height = 7, units = "in", res = 600
  )
  print(p_dot)
  dev.off()

  # ---- gene ratio plot ----
  plot_df <- enrich_df %>%
    slice_head(n = n_show) %>%
    mutate(
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
    file.path(save_dir, paste0(prefix, "_gene_ratio.tiff")),
    width = 9, height = 7, units = "in", res = 600
  )
  print(p_ratio)
  dev.off()

  png(
    file.path(save_dir, paste0(prefix, "_gene_ratio.png")),
    width = 9, height = 7, units = "in", res = 600
  )
  print(p_ratio)
  dev.off()
}
