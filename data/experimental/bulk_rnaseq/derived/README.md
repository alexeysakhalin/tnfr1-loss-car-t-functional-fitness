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
this snapshot before rendering the figures.
