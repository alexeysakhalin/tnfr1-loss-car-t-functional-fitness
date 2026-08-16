# Frozen analysis resources

`CAR_T_state_signatures.csv` contains ten fixed C0-C9 expression signatures,
with 20 positively ranked genes per computational cluster. The membership was
frozen before the revised public-cohort scoring and is checksum-locked by the
repository tests.

`CAR_T_state_signature_concordance_v1.csv` is a release-validation contract for
rerunning the exploratory marker ranking under R 4.4.3, Seurat 5.5.1 and Matrix
1.7.2. It records, for every C0-C9 cluster, the exact overlap and exact
current-only/frozen-only genes observed in that reviewed environment. The R/05
release guard accepts only that complete row-for-row result; it is not a
percentage-similarity tolerance. Any additional, missing or different gene
requires a new manual marker/label review and a deliberately versioned contract.
The contract documents reproducibility of the ranking and does not modify the
frozen clinical-scoring signatures.

These are top-ranked genes from an exploratory cell-level marker analysis,
not ten sets of biological-replicate differential-expression findings. No
adjusted-p-value filter is applied when the fixed 20-gene membership is read.
In the frozen table, 198 of 200 rows have adjusted p-values below 0.05; KLRC1
in C0 and CHI3L2 in C9 do not. This does not change membership, but it prevents
the signatures from being described as 200 statistically significant markers.

The single-cell experiment has one targeted count matrix per conditioning
label. Marker p-values and cluster fractions are therefore descriptive at the
cell/computational level and do not provide biological-replicate inference or
absolute viable-cell phenotype frequencies.

R/05 retains all QC-passing C0-C10 cells in the descriptive Figure 4A-B
denominator. C10 is reported neutrally as a small cytokine/IFN-responsive
cluster outside the historical C0-C9 frozen/transferred mapping; it is not
classified as a contaminant. Aggregate C10 marker, cluster-count and
filtering-QC tables are written to
`results/targeted-singlecell-diagnostics/` before the concordance guard. These
diagnostics remain available if the guard stops a release; no cell-level count
matrix is copied there.

`Figure_5F_curated_gene_sets.csv` defines the prespecified display categories
for the exploratory aggregate CheckMate gene-model plot. Identifier maps are
checksum-pinned snapshots used during cohort preparation.
