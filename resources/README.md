# Versioned analysis resources

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

`targeted_singlecell_cluster_annotations_v1.tsv` is the human-readable,
versioned annotation record for all tumor-co-culture C0-C10 and independently
clustered repeated-stimulation C0-C5 states. It records the earlier release
label, the property-based submission label, QC-passing cell count, defining
markers, marker-table source, interpretive literature and properties that were
not measured. Annotation was manual and post hoc after unsupervised Louvain
clustering; no reference atlas, label transfer or automated classifier was
used. Masopust et al. (2026,
[doi:10.1038/s41577-025-01238-2](https://doi.org/10.1038/s41577-025-01238-2))
is the nomenclature and reporting framework, not a source of transferred
cluster identities.

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

`IMvigor210_expression_semantic_contract_v1.json` contains no sample or feature
identifiers and no expression cells. It fixes the complete matrix dimensions,
ordered-identifier hashes, structural-zero count, direct-formula tolerance and
the unambiguous UTF-8/LF framing used by the streaming verifier. Fixed6
`Decimal`/`ROUND_HALF_UP` semantics are the required compatibility gate;
fixed7 and fixed8 digests are retained only to diagnose harmless rendering
differences. The canonical whole-file SHA-256 identifies one historical CSV
rendering for provenance and is not an expression acceptance criterion. The
same contract requires the preparer to replace each accepted cell by its fixed6
scaled value before mapping, aggregation, ranking or output. Compatibility is
therefore analysis-equivalent, including under adversarial sub-six-decimal
changes.

`IMvigor210_fixed6_canonicalization_impact_v1.json` records the aggregate,
non-sample-identifying historical A/B check. It names only two prespecified
genes used in the manuscript-level median comparison. All 61,944
selected-expression keys were identical;
59,954 displayed expression values moved by at most
`4.999993601373376e-7`, while zero rank percentiles changed. Figure 5C scores,
Supplementary Figure S6 score orderings/Spearman results and all tested
selection or significance thresholds were unchanged. A regenerated density
plot can differ at the file-byte level because it uses the rounded expression
coordinates, but its numerical interpretation is unchanged.
