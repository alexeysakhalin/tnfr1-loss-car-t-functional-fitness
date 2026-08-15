# Bulk RNA-seq analysis specification

## Input contract

The analysis starts from
`data/experimental/bulk_rnaseq/gene_counts.tsv.gz`. Rows are unique source
`Gene_ID` values and columns are the 24 samples declared in
`sample_metadata.tsv`. Counts are complete, non-negative integers. Source
annotations are stored separately so the count matrix remains machine-readable
and identifier-stable.

The only repeated `Gene_Symbol` is `TRNAV-CAC`, assigned by the source to
Gene_ID 107985614, 107985615 and 107985753. Before differential-expression
model fitting, counts for these three identifiers are summed. All other symbols
map one-to-one to Gene_ID. The canonical Gene_ID table is retained unchanged;
the aggregation is an explicit analysis transformation.

## Experimental factors

The matrix contains WT and TNFR1-KO1 HeLa CD19-positive cells. The source name
`T6` is mapped to TNFR1-KO1. Treatments were 50 ng/mL TNF, 50 ng/mL IFNγ,
their combination, or matched untreated control for 48 hours. Each group has
three independent experimental replicates. Pairing of corresponding WT and
KO1 replicate numbers by experimental batch was not confirmed, so the primary
models are unpaired.

## Primary manuscript-reproduction models

The workflow uses PyDESeq2 0.5.4 with median-of-ratios size factors, the
parametric dispersion trend, Wald tests, Cook's filtering and independent
filtering. No LFC shrinkage is applied.

Each comparison is fitted as a separate two-group, six-sample model. This
retains the pairwise normalization and dispersion strategy used for the legacy
tables.

| Output stem | Model | Explicit contrast | Meaning of positive log2 fold change |
|---|---|---|---|
| `WT_IFNG_vs_control` | `~treatment` | IFNγ / control in WT | Higher after IFNγ |
| `WT_TNF_vs_control` | `~treatment` | TNF / control in WT | Higher after TNF |
| `WT_TNF_IFNG_vs_control` | `~treatment` | TNF+IFNγ / control in WT | Higher after TNF+IFNγ |
| `TNFR1_KO1_vs_WT_control` | `~genotype` | KO1 / WT in control | Higher in KO1 |
| `TNFR1_KO1_vs_WT_IFNG` | `~genotype` | KO1 / WT after IFNγ | Higher in KO1 |
| `TNFR1_KO1_vs_WT_TNF` | `~genotype` | KO1 / WT after TNF | Higher in KO1 |
| `TNFR1_KO1_vs_WT_TNF_IFNG` | `~genotype` | KO1 / WT after TNF+IFNγ | Higher in KO1 |

## Prespecified interaction analysis

With `--include-interaction`, all 24 samples are fitted with:

```text
~genotype + treatment + genotype:treatment
```

For each cytokine condition, the numeric contrast is:

```text
(KO1 - WT under treatment) - (KO1 - WT under control)
```

It therefore tests whether the treatment response differs by genotype. These
difference-in-differences results are not interchangeable with the
within-treatment KO1-versus-WT contrasts used in the figures. The pinned
24-sample workflow was run for the release snapshot and the three complete
interaction tables are retained under `data/experimental/bulk_rnaseq/derived/`.

## Output contract

Every result filename ends in `.unfiltered.tsv.gz`. Each file contains the same
46,425 sorted symbols. The workflow applies no export filter based on base mean,
fold change, p-value or adjusted p-value. Features that are all zero within a
particular six-sample subset are restored after model fitting with explicit
`NA` statistics.

| Column | Definition |
|---|---|
| `Gene_Symbol` | Unique analysis symbol |
| `source_gene_ids` | Semicolon-delimited source Gene_ID membership |
| `n_source_gene_ids` | Number of source Gene_ID rows summed |
| `total_count_in_model_subset` | Raw count sum over samples in the fitted model |
| `analysis_status` | `modelled`, `all_zero_in_subset`, or `not_returned_by_engine` |
| `baseMean` | PyDESeq2 mean of normalized counts for the fitted six- or 24-sample model |
| `log2FoldChange` | Numerator relative to denominator as declared above |
| `lfcSE` | Standard error of the log2 fold change |
| `stat` | Wald statistic; validated against `log2FoldChange / lfcSE` |
| `pvalue` | Unadjusted Wald-test p-value |
| `padj` | Benjamini-Hochberg adjusted p-value after independent filtering |

Expected primary output files are:

```text
results/bulk_rnaseq/WT_IFNG_vs_control.unfiltered.tsv.gz
results/bulk_rnaseq/WT_TNF_vs_control.unfiltered.tsv.gz
results/bulk_rnaseq/WT_TNF_IFNG_vs_control.unfiltered.tsv.gz
results/bulk_rnaseq/TNFR1_KO1_vs_WT_control.unfiltered.tsv.gz
results/bulk_rnaseq/TNFR1_KO1_vs_WT_IFNG.unfiltered.tsv.gz
results/bulk_rnaseq/TNFR1_KO1_vs_WT_TNF.unfiltered.tsv.gz
results/bulk_rnaseq/TNFR1_KO1_vs_WT_TNF_IFNG.unfiltered.tsv.gz
```

The optional interaction outputs use the prefix
`interaction_TNFR1_KO1_vs_WT_`. `analysis_manifest.tsv` records model,
contrast direction, row count, status counts and output checksum.
`run_metadata.json` records input hashes and direct runtime versions, and
`environment.freeze.txt` records the complete installed Python environment.

The `results/bulk_rnaseq/figure_inputs/` adapters contain snake-case schemas
for the figure scripts. Both use the unambiguous column `gene_symbol`; no gene
symbol is labelled as a Gene_ID. WT effect columns explicitly state
`treatment_vs_untreated`, and KO effect columns explicitly state
`ko1_vs_wt`. The combined WT table is supplied both without export filtering
and as an explicitly named `fdr05` view. The latter is a convenience view and
is not an acceptable input for reconstructing a volcano plot's non-significant
background.

The manuscript figure scripts read the validated adapters committed under
`data/experimental/bulk_rnaseq/derived/`. A count-level rebuild first writes
new adapters under `results/bulk_rnaseq/figure_inputs/`; the release verifier
must accept that rebuild before the committed adapters are used to render the
figures. The R scripts use the snake-case effect columns declared above; the
legacy workbooks and their historical column aliases are not mixed with these
adapters.

`requirements-bulk-rnaseq.txt` records the complete exact Python 3.12.13
environment used for the release run. The package snapshot is retained as
`environment.freeze.txt` with the derived results. Each CI rebuild also writes
`runtime_provenance.json`, which records the operating platform, Python and
zlib versions, NumPy build configuration, and numerical-library thread
environment. Runtime provenance is diagnostic and is not required to match a
different host.

## Release-outcome check

`scripts/verify_bulk_rnaseq_release.py` keeps schema, semantic keys and row
order, text and integer fields, `NA` masks, analysis metadata, and manifest
provenance exact. The generated manifest hashes must match the generated files.
For every release table, the base-mean and fold-change threshold masks and the
combined DEG categories (`baseMean >= 30`, adjusted p-value `< 0.05`, and
log2 fold change `> 1` or `< -1`) must be identical. Figure 1C and
Supplementary Figure S2D Venn membership and the significance/direction calls
for ICAM1, MLKL, GSDME, and IRF1 in the interaction analyses must also be
identical. Both the regenerated and committed tables are independently checked
for the Wald identity `stat = log2FoldChange / lfcSE`. For the 21 genes labelled
in the volcano plots, regenerated log2 fold changes must be within an absolute
difference of `0.001` of the committed value wherever estimable.

Small changes in numerically unstable model tails can occur across BLAS and
libm implementations even with the locked package environment. Raw numeric
deltas and standalone adjusted-p-value threshold flips are therefore recorded
in `release_comparison_report.json`; they do not fail the rebuild unless one of
the exact scientific-outcome contracts above changes.
