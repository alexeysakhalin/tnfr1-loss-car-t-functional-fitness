# Figure reproduction map

This document maps each manuscript panel to its executable script, required
input and expected output. Supplementary numbering assumes the acute
ICAM1/conjugate experiment is Supplementary Figure S3.

## Experimental transcriptomic panels

| Panel | Script | Version-controlled input | Expected output |
|---|---|---|---|
| Figure 1B | `scripts/run_bulk_rnaseq_pydeseq2.py`; `R/03_figure_1_B_C.R` | `data/experimental/bulk_rnaseq/`; complete WT cytokine adapters | `figures/figure_1/Figure_1B_triptych.{png,tiff}` |
| Figure 1C | same workflow | complete WT cytokine adapters | `figures/figure_1/Figure_1C_upregulated_overlap.{png,tiff}` |
| Figure 2B | `scripts/run_bulk_rnaseq_pydeseq2.py`; `R/04_figure_2_B_suppl_S2D.R` | complete within-treatment TNFR1-KO1-versus-WT adapters | `figures/figure_2/Figure_2B_triptych.{png,tiff}` |
| Figure 2C | same workflow | ICAM1 and IRF1 rows from the complete within-treatment adapter | `figures/figure_2/Figure_2C_ICAM1_IRF1_effects.{png,tiff}`; exact estimates in `results/figure_2/Figure_2C_ICAM1_IRF1_effects.tsv` |
| Supplementary Figure S2D | same workflow | combined four-stratum adapter; Venn uses TNF, IFN-gamma and TNF+IFN-gamma | `figures/figure_2/Supplementary_Figure_S2D_downregulated_overlap.{png,tiff}` |
| Figure 4A-B | `R/05_figure_4_AB_suppl_S5A.R` | `WT_targeted_counts.tsv.gz`, `KO1_targeted_counts.tsv.gz`, `KO2_targeted_counts.tsv.gz` | `figures/Figure_4A_UMAP_clusters_annotated.png`; `figures/Figure_4B_cluster_fraction_facet.png` |
| Supplementary Figure S5A | `R/05_figure_4_AB_suppl_S5A.R` | same WT/KO matrices | `figures/Figure_4A_UMAP_4panel.{png,tiff}` |
| Figure 5A | `R/05_figure_4_AB_suppl_S5A.R` | `TCR_targeted_counts.tsv.gz`, 5,662 CD3-positive cells | `figures/Figure_5A_TCR_UMAP_clusters.png`; `figures/Figure_5A_TCR_cluster_composition.png` |

The bulk panels use complete regenerated result universes. Figure 1B/1C DEG
sets require baseMean >=30, adjusted P <0.05 and the stated fold-change
threshold. Figure 2B uses within-treatment TNFR1-KO1-versus-WT contrasts, not
the genotype-by-treatment interaction coefficients. Figure 2C shows the same
within-treatment model estimates for ICAM1 and IRF1 with unadjusted 95% Wald
confidence intervals and within-contrast BH-adjusted P values; it contains no
between-treatment tests. Source label `T6` maps to manuscript clone
`TNFR1-KO1`.

For the targeted single-cell panels, archive before/after QC counts, nonempty
marker tables and `sessionInfo()`. Manual cluster labels must be checked against
the exported marker profiles because numerical cluster identifiers may permute
between runs.

## Clinical-context panels

| Panel | Script | Local input | Expected output |
|---|---|---|---|
| Figure 5C | `R/06_figure_5_C.R` | prepared four-cohort expression and metadata tables | `figures/Figure_5C_bulk_Tcell_expression_score_rank_based.png` |
| Figure 5D | `R/07_figure_5D.R` | nivolumab-only CheckMate tables | `figures/Figure_5D_KM_ccRCC_Tcell_score.png` |
| Figure 5E | `R/08_figure_5E_S6D_S6E.R` | nivolumab-only CheckMate tables; frozen C0-C9 signatures | `figures/Figure_5E_KM_ccRCC_C6.png` |
| Figure 5F | `R/09_figure_5F_S6G.R` | tracked aggregate global gene models | `results/figure5F/figures/Figure_5F_C6_adjusted_gene_level_volcano.png` |
| Figure 5G | `R/10_figure_5G.R` | nivolumab-only CheckMate tables; frozen C0-C9 signatures | `figures/Figure_5G_B2M_ICAM1_signature_associations.png` |
| Supplementary Figure S6A,C | `R/02_supplementary_S6A_S6C.R` | prepared four-cohort expression and metadata tables | `figures/Supplementary_Figure_S6A_*`; `figures/Supplementary_Figure_S6C_*` |
| Supplementary Figure S6B,F | `R/12_supplementary_S6B_S6F.R` | prepared four-cohort expression and metadata tables | `results/supplementary_S6/Supplementary_Figure_S6B_*`; `Supplementary_Figure_S6F_*` |
| Supplementary Figure S6D,E | `R/08_figure_5E_S6D_S6E.R` | nivolumab-only CheckMate tables | `figures/Supplementary_Figure_S6D_*`; `figures/Supplementary_Figure_S6E_*` |
| Supplementary Figure S6G | `R/09_figure_5F_S6G.R` | tracked aggregate group-balance table | `results/figure5F/figures/Supplementary_Figure_S6G_C6_group_balance.png` |

All CheckMate clinical-context panels use the 181 nivolumab-treated tumors.
The Figure 5 and Supplementary Figure S6 composites must be assembled entirely
from these outputs; pooled-arm panels are not compatible with this analysis.

## DepMap and flow-cytometry panels

| Panel | Scope | Reproduction requirement |
|---|---|---|
| Supplementary Figure S1B | repository workflow | `R/11_supplementary_1B.R` renders the tracked 1,591-model derivative; exact contracts are in `docs/DEPMAP_S1B.md` |
| Supplementary Figure S3 | separate flow-cytometry workflow | raw FCS files, compensation/gating workspace, live-cell denominator, E:T ratio, timing, replicate/donor map and prespecified WT-versus-KO contrasts |

The DepMap expression archive and `Model.csv` are a confirmed same-release
DepMap Public 25Q2 pair downloaded from the portal's **All Data** page. The
checksum-pinned identities and release assertion are recorded in the
provenance file; no release-specific Figshare DOI is assigned to 25Q2.

## Final assembly checks

1. Run the scripts from a clean clone and archive `sessionInfo()`, console logs
   and source checksums.
2. Use the 600-dpi TIFF files or journal-approved vector exports for the final
   article; retain PNG files for visual checking.
3. Verify every reported number against the generated TSV/CSV result tables.
4. Do not combine panels generated from different cohort definitions or source
   versions in one composite.
