# Figure reproduction and manuscript replacement map

This map follows the supplementary numbering after insertion of the acute
ICAM1/conjugate experiment as Supplementary Figure S3. It distinguishes panels
that can be regenerated now from panels that remain locked by missing source
data or metadata.

## Experimental panels

| Manuscript panel | Script and output | Input | Manuscript action |
|---|---|---|---|
| Figure 1B | `scripts/run_bulk_rnaseq_pydeseq2.py`, then `R/03_figure_1_B_C.R`; final composite: `figures/figure_1/Figure_1B_triptych.png` and `Figure_1B_triptych.tiff`; audit panels: `Figure_1B_TNF_IFNg_volcano.*`, `Figure_1B_IFNg_volcano.*`, `Figure_1B_TNF_volcano.*` | Version-controlled count matrix, sample metadata and validated complete adapter under `data/experimental/bulk_rnaseq/derived/` | Use the 600-dpi LZW TIFF triptych in the manuscript, in the fixed order TNF+IFNγ, IFNγ, TNF. Every modelled gene with a finite adjusted p-value is shown; the common symmetric x range is rounded outward from all plotted effects. Adjusted P values at or below 1×10^-300 are displayed at -log10(adjusted P)=300, as stated on the shared y-axis. DEG colour uses baseMean ≥30, absolute log2 fold change >1 and adjusted P <0.05. Labels are limited a priori to ICAM1, IRF1, AIM2, CASP1, CASP4, MLKL, BAK1, CASP7, FAS and CASP8. Positive log2 fold change means higher expression after treatment. |
| Figure 1C | Same workflow; `figures/figure_1/Figure_1C_upregulated_overlap.*` | Same complete regenerated adapters; the overlap sets use baseMean ≥ 30, adjusted P < 0.05 and log2 fold change > 1 | Replace together with Figure 1B and verify the exported set/intersection counts. |
| Figure 2B | `scripts/run_bulk_rnaseq_pydeseq2.py`, then `R/04_figure_2_B_suppl_S2D.R`; final composite: `figures/figure_2/Figure_2B_triptych.png` and `Figure_2B_triptych.tiff`; audit panels: `Figure_2B_TNF_IFNg.*`, `Figure_2B_IFNg.*`, `Figure_2B_TNF.*` | Complete regenerated TNFR1-KO1-versus-WT adapters; T6 is author-confirmed as TNFR1-KO1 | Use the 600-dpi LZW TIFF triptych in the manuscript, in the fixed order TNF+IFNγ, IFNγ, TNF. The common symmetric x range is rounded outward from all plotted effects. Adjusted P values at or below 1×10^-300 are displayed at -log10(adjusted P)=300, as stated on the shared y-axis. Labels are limited a priori to ICAM1, MLKL, GSDME and IRF1. Positive log2 fold change means higher expression in TNFR1-KO1 within the indicated condition. |
| Supplementary Figure S2D | Same workflow; `figures/figure_2/Supplementary_Figure_S2D_downregulated_overlap.*` | The combined adapter contains all four strata for universe/QC checks; the Venn analysis uses the TNF, IFNγ and TNF+IFNγ contrasts | Replace together with Figure 2B and verify the exported set/intersection counts. |
| Figure 4A-B | `R/05_figure_4_AB_suppl_S5A.R`; `figures/Figure_4A_UMAP_clusters_annotated.png`, `Figure_4B_cluster_fraction_facet.png` | WT/KO1/KO2 targeted counts from CD3+ sheets | Replace both panels after a clean Seurat run. Cluster fractions describe computationally retained cells, not absolute live-cell phenotype frequencies. |
| Supplementary Figure S5A | `R/05_figure_4_AB_suppl_S5A.R`; `figures/Supplementary_Figure_S5A_UMAP_by_sample.png` | Same WT/KO inputs | Replace after the same clean Seurat run. |
| Figure 5A | `R/05_figure_4_AB_suppl_S5A.R`; `figures/Figure_5A_TCR_UMAP_clusters.png`, `Figure_5A_TCR_cluster_composition.png` | TCR targeted counts from the CD3+ sheet, n=5,662 | Replace after a clean Seurat run. The six labels are independent of the tumor-co-culture C0-C9 labels. |

The source TCR workbook's `all cells` sheet contains three additional
CD3-negative cells and is not used. All four single-cell matrices contain 259
targeted genes; they must be described as targeted single-cell mRNA profiling.
The clean run must produce before/after QC counts and nonempty marker tables,
and the authors must review marker profiles against every manual label.
Reappearance of cluster IDs 0-N is not sufficient because graph-cluster
numbers can permute between runs.

The original integer count matrix, sample manifest and count-to-DE code are
version-controlled. `scripts/run_bulk_rnaseq_pydeseq2.py` generates three WT
cytokine-versus-control contrasts, four within-condition
TNFR1-KO1-versus-WT contrasts and, when requested, three separate formal
genotype-by-treatment interaction contrasts. `R/03` and `R/04` accept only the
complete regenerated adapters; the interaction results do not replace the
within-condition Figure 2 contrasts.

The numerical display contracts are written to
`results/figure_1/Figure_1B_volcano_output_contract.tsv` and
`results/figure_2/Figure_2B_volcano_output_contract.tsv`. They record the fixed
panel order, plottable row counts, observed effect extrema, common x limits,
P-value floor, y-cap counts and the labels present in each panel. The figure
workflow checks these tables against the complete adapters and verifies the
pixel dimensions, 600-dpi metadata and LZW compression of every PNG/TIFF pair.

## Clinical-context panels

| Manuscript panel | Script and output | Required local input | Manuscript action |
|---|---|---|---|
| Figure 5C | `R/06_figure_5_C.R`; `figures/Figure_5C_bulk_Tcell_expression_score_rank_based.png` | Locally prepared four-cohort tables | Replace the current temporary panel. The y-axis is a rank-based bulk T-cell expression score, not infiltration or a cell fraction. |
| Figure 5D | `R/07_figure_5D.R`; `figures/Figure_5D_KM_ccRCC_Tcell_score.png` | Local nivolumab-only CheckMate table | Current numerical panel is consistent with the reference results; regenerate in the final clean run for uniform styling. |
| Figure 5E | `R/08_figure_5E_S6D_S6E.R`; `figures/Figure_5E_KM_ccRCC_C6.png` | Same n=181 CheckMate table and frozen C0-C9 signatures | Replace. Report Cox and log-rank p-values separately; neither multiplicity-adjusted split result is significant. |
| Figure 5F | `R/09_figure_5F_S6G.R`; `results/figure5F/figures/Figure_5F_C6_adjusted_gene_level_volcano.png` | Version-controlled aggregate global gene models | Replace; the current manuscript panel predates the frozen aggregate model table. |
| Figure 5G | `R/10_figure_5G.R`; `figures/Figure_5G_B2M_ICAM1_signature_associations.png` | n=181 CheckMate table and frozen C0-C9 signatures | Replace. This is an exploratory 30-model family with global BH correction. |
| Supplementary Figure S6A,C | `R/02_supplementary_S6A_S6C.R`; `figures/Supplementary_Figure_S6A_*`, `Supplementary_Figure_S6C_*` | Locally prepared four-cohort tables | Regenerate all panels. S6C is descriptive and its visually selected TNFRSF1A cutoff is not an exclusion rule. |
| Supplementary Figure S6B,F | `R/12_supplementary_S6B_S6F.R`; `results/supplementary_S6/Supplementary_Figure_S6B_*`, `Supplementary_Figure_S6F_*` | Locally prepared four-cohort tables | Replace. Twelve Spearman tests use one BH family; trend lines must not be interpreted as the Spearman model. |
| Supplementary Figure S6D,E | `R/08_figure_5E_S6D_S6E.R`; `figures/Supplementary_Figure_S6D_*`, `Supplementary_Figure_S6E_*` | n=181 CheckMate table | Replace because the current embedded values are stale. |
| Supplementary Figure S6G | `R/09_figure_5F_S6G.R`; `results/figure5F/figures/Supplementary_Figure_S6G_C6_group_balance.png` | Version-controlled aggregate balance table | Replace; correct group sizes are Low n=91 and High n=90. |

For the clean article build, replace the complete Figure 5 C-G block and the
complete Supplementary Figure S6 A-G composite. This avoids mixing panels made
from the old pooled n=311 workflow with panels from the corrected n=181
nivolumab-only workflow.

## Source-locked panels

| Panel | Status | Required resolution |
|---|---|---|
| Supplementary Figure S1B | Clean-clone rendering and numerical contract complete; metadata-release wording pending | `scripts/prepare_depmap_s1b.py` validates the raw expression ZIP and `Model.csv`, selects one default profile per `ModelID`, and writes a deterministic tracked derivative containing 1,591 eligible cell-line models. `R/11_supplementary_1B.R` renders one manuscript-ready scatter directly from that derivative and rechecks all flags, quadrants and counts. The fixed results are 1,003/1,591 RIPK3 below threshold, 1,172/1,591 NLRP3 below threshold and 749/1,591 both below threshold. `TissueOrigin` is empty and is not used. The expression source is fixed as DepMap Public 25Q2; the `Model.csv` release identity and release-pair status remain machine-locked as unverified and cannot be overridden at render time. Before submission, verify the metadata release or disclose that it was separately checksum-pinned with unverified release identity. Use “eligible DepMap cell-line models,” not “human cell lines,” and do not present the percentages as prevalence or pathway competence. |
| Supplementary Figure S3 | Not generated by this repository | Supply raw FCS files, compensation/gating workspace, live/dead gating and denominator, E:T ratio, timing, cytokine-exposure timing, replicate/donor map and exact WT-versus-KO contrasts. Replace `Delta MFI normalized fold change` with a defined recovered dual-positive-event metric. |

The acute S3 assay supports the statement that no reduction was detected in
the recovered dual-positive readout for TNFR1-KO targets under the tested
conditions. It does not establish immune-synapse architecture, contact
duration, serial engagement or LFA-1 dependence.

## Final assembly checks

1. Run the scripts from a clean clone with the release `renv.lock`.
2. Archive `sessionInfo()`, console logs and source-table checksums.
3. Use the 600-dpi TIFF files or journal-approved vector exports in the
   manuscript; PNG files are for visual comparison.
4. Compare every number in the figure, legend, Results, Discussion, Table 1
   and rebuttal against the generated TSV/CSV result table.
5. Do not mix old and regenerated subpanels in one composite.
6. Do not restore the removed legacy reference PNGs; add new checksum-pinned
   references only after the matching QC, marker, count and `sessionInfo()`
   records have been archived.
