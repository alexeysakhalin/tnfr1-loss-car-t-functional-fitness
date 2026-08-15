# Validated bulk RNA-seq release tables

These files are the complete, unfiltered outputs used by the manuscript bulk
RNA-seq figures and the prespecified genotype-by-treatment analysis. They were
generated from the version-controlled Gene_ID-level integer count matrix with
Python 3.12.13 and PyDESeq2 0.5.4 using two CPUs.

The two `figure_*` adapters contain the seven primary six-sample contrasts in
the exact schemas consumed by `R/03_figure_1_B_C.R` and
`R/04_figure_2_B_suppl_S2D.R`. The three `interaction_*` files contain the
formal difference-in-differences contrasts

`(TNFR1-KO1 treated - WT treated) - (TNFR1-KO1 control - WT control)`.

Every contrast retains all 46,425 gene symbols. All-zero features have
`baseMean=0` and missing inferential fields; no significance or expression
filter was applied before export. `analysis_manifest.tsv`,
`run_metadata.json`, `environment.freeze.txt` and `SHA256SUMS` record the
validated release run. The automated rebuild compares regenerated outputs to
this snapshot before rendering the figures. Semantic identifiers, annotations,
integer counts, missing-value locations, table order, run metadata, the complete
environment freeze, and all non-hash manifest fields must match exactly. Output
SHA-256 values are validated against each regenerated file and are excluded
only from the cross-platform manifest comparison.

Raw floating-point differences are recorded as diagnostics because
numerically unstable model tails can vary across BLAS/libm implementations.
They are not subject to a global closeness threshold. Instead, the verifier
requires exact membership for `baseMean >= 30`, log2 fold change `> 1` and
`< -1`, the combined up/down DEG categories at adjusted p-value `< 0.05`, and
the Figure 1C and Supplementary Figure S2D Venn sets. Significance and direction
must remain exact for ICAM1, MLKL, GSDME, and IRF1 in every interaction table.
For the 21 genes labelled in the manuscript volcano plots, regenerated log2
fold changes must additionally be within an absolute difference of `0.001` of
the release value wherever estimable. Every regenerated and release table is
also checked internally for `Wald = log2 fold change / standard error`.

Standalone adjusted-p-value threshold flips are reported and are accepted only
when they leave the combined DEG category unchanged. CI preserves the JSON
comparison report, regenerated results, and runtime provenance even when
verification fails. Figures are rendered from the committed release adapters
only after the outcome gate passes.
