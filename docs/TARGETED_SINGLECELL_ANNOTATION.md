# Targeted single-cell cluster annotation and nomenclature

## Scope

The study used targeted single-cell mRNA profiling with the 259-gene BD
Rhapsody Human T-cell panel. The labels therefore describe expression
properties observed within that panel; they are not whole-transcriptome cell
types.

The tumor-co-culture and repeated-CD3/CD28-stimulation datasets were analyzed
independently. A cluster identifier such as C0 has meaning only within its own
dataset. Identically numbered clusters are not assumed to be homologous.

## Annotation method

Clusters were generated de novo with the workflow in
`R/05_figure_4_AB_suppl_S5A.R`: LogNormalize with a scale factor of 10,000,
regression of `nCount_RNA`, 30 principal components, Louvain clustering on
dimensions 1-20 at resolution 0.6, and UMAP on dimensions 1-20 with seed 12,345.

Positive cluster-enriched genes were ranked with
`Seurat::FindAllMarkers` using the Wilcoxon rank-sum test,
`only.pos = TRUE`, `logfc.threshold = 0.25`, and `min.pct = 0.10`.
Labels were assigned manually after clustering from these marker profiles. No
reference-atlas label transfer, SingleR, Azimuth, CellTypist, or other
automated cell-type classifier was used.

The computational framework is Seurat v5
([Hao et al., 2024](https://doi.org/10.1038/s41587-023-01767-y)).
The naming policy follows the property-based reporting recommendations of
[Masopust et al., 2026](https://doi.org/10.1038/s41577-025-01238-2).
Masopust et al. is a nomenclature and reporting framework, not a reference
atlas and not the source of transferred labels.

Terms such as *associated*, *high*, and *-like* denote marker-defined
transcriptional resemblance. They do not establish lineage, ontogeny, antigen
specificity, protein function, cytotoxic capacity, tissue residence, memory
potential, or exhaustion. Cell-level marker P values are descriptive and are
not treated as biological-replicate inference.

## Versioned annotation table

The complete record is
[`resources/targeted_singlecell_cluster_annotations_v1.tsv`](../resources/targeted_singlecell_cluster_annotations_v1.tsv).
Its `analysis_source_commit` field anchors the unchanged clustering, cell
membership and count-matrix release used for annotation review; the Git tag or
commit containing this manifest separately identifies the nomenclature and
projection release.
For each tumor-co-culture C0-C10 and repeated-stimulation C0-C5 cluster it
records:

- QC-passing cell count;
- the previous release label and the property-based submission label;
- defining marker anchors;
- marker-table source;
- annotation method and assay scope;
- the role of interpretive literature;
- properties that were not measured;
- whether the cluster participates in a transferred signature family.

The supporting marker outputs are
`Supplementary_Table_S5_top_markers_per_cluster.xlsx`,
`Supplementary_Table_S5_C10_markers.tsv`, and
`Supplementary_Table_TCR_top_markers_per_cluster.xlsx` in the
`targeted-singlecell-r05-results` workflow artifact.

### Annotation summary

The table below is a compact reader-facing view of the versioned manifest.
Marker anchors are descriptive evidence from the targeted panel, not a claim
that every cell in the cluster expresses every listed gene.

#### Tumor-co-culture dataset

| Cluster | QC cells | Submission label | Defining marker anchors |
|---|---:|---|---|
| C0 | 1,542 | KLRB1/LGALS3-associated activated T-cell state | LGALS3, KLRB1, CD69, CXCR4, ITGAE |
| C1 | 1,448 | CD4/LAG3-associated activated T-cell state | CD4, LAG3, TNFRSF4, TNFRSF18, STAT3 |
| C2 | 903 | CXCR5/IL13/CCR4-associated activated T-cell state | CXCR5, IL13, CCR4, CCR8, CCL1, CSF2 |
| C3 | 891 | CXCR6-associated cytotoxic activated state | CXCR6, CCL3, GZMB, PRF1, CCR5 |
| C4 | 812 | CD8/ZNF683-associated cytotoxic state | ZNF683, CD8A, CD8B, GZMH, NKG7, KLRK1 |
| C5 | 594 | TRDC-high γδ-associated cytotoxic state | TRDC, TARP-REFSEQ, ZBTB16, NKG7, NCR3 |
| C6 | 516 | CXCL13-associated cycling T-cell state | CXCL13, TK1, MKI67, AURKB, TOP2A, UBE2C |
| C7 | 453 | Cycling effector-gene-high T-cell state | MKI67, TOP2A, AURKB, UBE2C, CCNB1, HMGB2 |
| C8 | 396 | IL9-high activated T-cell state | IL9, TNFRSF8, IL5, STAT6, CCR4, CCR8 |
| C9 | 133 | TCF7/IL7R-high early-memory-associated state | TCF7, IL7R, LEF1, CCR7, SELL, CD27 |
| C10 | 33 | Small cytokine/IFN-response-high cluster | CXCL9, CXCL10, IL6, IL15, IFNGR1, OAS1, STAT1 |

#### Repeated CD3/CD28-stimulation dataset

| Cluster | QC cells | Submission label | Defining marker anchors |
|---|---:|---|---|
| C0 | 1,298 | Mixed CD4/KLRB1-associated activated state | FOXP3, KLRB1, CD4, CTLA4, HAVCR2, PRF1 |
| C1 | 1,034 | Cycling T-cell state I | HMMR, CCNB1, UBE2C, MKI67, AURKB, TOP2A |
| C2 | 1,015 | Cytokine-expressing effector state | CCL1, CCL4, CSF2, IFNG, XCL1, GZMB, FASLG |
| C3 | 833 | CD8/TRDC-associated cytotoxic state | CD8A, CD8B, ZNF683, TRDC, KLRK1, NKG7, CTSW |
| C4 | 806 | Cycling T-cell state II | HMMR, UBE2C, TOP2A, AURKB, MKI67, CCNB1 |
| C5 | 444 | CCR7/IL7R/HLA-II-associated state | CCR7, IL7R, LEF1, HLA-DRA, HLA-DQA1, HLA-DPA1 |

## CXCL13-associated tumor-co-culture C6

Tumor-co-culture C6 remains the **CXCL13-associated cycling T-cell state**.
Its marker profile includes CXCL13 together with TK1, MKI67, AURKB, TOP2A,
UBE2C, HMGB2, TYMS, HMMR, and PTTG2. The term *CXCL13-associated* is deliberate:
not every C6 cell has detectable CXCL13, and the label does not establish a
stable lineage or exhaustion program. Literature on CXCL13-positive
tumor-reactive T-cell states provides interpretive context
([Liu et al., 2022](https://doi.org/10.1038/s43018-022-00433-7)); it was not
used as an annotation classifier.

## Exploratory projection into the repeated-stimulation dataset

The frozen tumor-co-culture C6 signature contains 20 genes. Before inspecting
projection results it was partitioned into:

- cycle-associated: TK1, MKI67, AURKB, TOP2A, UBE2C, HMGB2, TYMS, HMMR, PTTG2;
- non-cycle/context: CHI3L2, CXCL13, FOXP1, CD70, IER5, IL23R, JUN, CXCR4, FAS,
  CCR7, CD4.

All 20 genes are present in the shared targeted panel. For each cell and each
gene set, genes are ranked within the fixed 259-gene panel with the
highest-expressed gene assigned rank 1 and tied values assigned their average
rank. The panel-relative Mann-Whitney rank-AUC score is

`score = 1 - U / [m(259 - m)]`,

where `m` is the number of signature genes and `U` is their rank-sum
statistic. Full, cycle, and non-cycle scores are reported separately because a
high full score can be driven by proliferation genes alone.

The projection is descriptive. It does not alter the repeated-stimulation
clustering or labels and does not identify repeated-stimulation C1 or C4 as
tumor-co-culture C6 cells. The public artifact contains only the rendered
figure, six-row cluster summary, and 20-row gene-coverage table; it excludes
cell barcodes, cell-level scores, expression matrices, and the Seurat object.
