# Experimental analysis tables

This directory contains compact, non-identifying, author-generated tables used
by the experimental figure scripts. The eight legacy result/count workbooks
are not used directly by R and are not committed here. Their canonical names, byte
sizes and SHA-256 values are recorded in `experimental_data_manifest.tsv`.
The manifest uses stable canonical deposition names. Timestamp suffixes added
by file transfer are removed from those names; the byte-level source identity
is preserved by SHA-256. The longer supplied IFN filename containing
`no_zero_pvals` is normalized because two p-values are zero from numerical
underflow. Deposit the unchanged source workbooks under the manifest's stable
canonical names; SHA-256, rather than a transient upload suffix, defines their
identity.

To rebuild these tables from the source workbooks, pass either a flat directory
or the root of the deposited archive containing `bulk_rnaseq/` and
`targeted_single_cell/` subdirectories:

```bash
python scripts/prepare_experimental_analysis_tables.py --input-dir /path/to/workbooks
```

The preparation script validates workbook structure, sample balance, gene and
cell identifiers, numeric ranges, duplicate records and count-matrix
integrity. Compressed outputs are written deterministically (`gzip` timestamp
zero) and replaced atomically.

## Cytokine-versus-untreated differential expression

`hela_cytokine_significant_differential_expression.tsv.gz` combines the TNF,
IFN-gamma and TNF+IFN-gamma source workbooks. Every supplied row has an
adjusted p-value below 0.05, so this is explicitly a significant-gene table,
not a full DESeq2 result table.

The source workbooks contain a contrast-orientation inconsistency: the sign of
`stat` is opposite to `log2FoldChange` in every row, although its magnitude and
the two-sided p-values are internally consistent. The canonical table retains
the treatment-versus-untreated fold-change orientation and recalculates the
Wald statistic as `log2FoldChange / lfcSE`. Source workbooks remain unchanged.

The significant-gene table is retained as a documented legacy input and QC
cross-check. It is not used to draw the final Figure 1B or Figure 1C. The
validated Macrogen integer count matrix, sample metadata, complete PyDESeq2
workflow and full unfiltered release tables are provided under `bulk_rnaseq/`;
R03 reads the version-controlled complete adapter from that directory.

## Matched T6-versus-WT RNA-seq

The four `hela_t6_vs_wt_*_differential_expression.tsv.gz` files contain the
condition-specific legacy tables from the matched workbook. The associated
design contains 24 samples: HeLa and T6, four conditions, three independent
experiments per cell-line/condition combination. They are retained for source
cross-checks; final Figure 2B and Supplementary Figure S2D use the count-level
rerun under `bulk_rnaseq/derived/`.

`T6` is the internal cell-line identifier encoded in the workbook. The authors
confirmed on 2026-08-14 that T6 corresponds to manuscript clone TNFR1-KO1.
The canonical table fields preserve `T6` for source-level provenance, while
figure labels use `TNFR1-KO1`. The mapping is recorded in
`experimental_sample_aliases.tsv`.

## Targeted single-cell count matrices

The four files under `singlecell/` are integer count matrices from the `CD3+
cells` sheet of each source workbook:

| Sample | Cells | Genes |
|---|---:|---:|
| WT | 3,743 | 259 |
| KO1 | 2,475 | 259 |
| KO2 | 1,831 | 259 |
| TCR repeated-stimulation dataset | 5,662 | 259 |

The TCR workbook's `all cells` sheet contains three additional CD3-negative
cells. They are intentionally excluded so that the stated T-cell analysis and
the executable input use the same population. These are targeted BD Rhapsody
panel data, not whole-transcriptome single-cell RNA sequencing.

## Redistribution

The canonical tables in this directory are author-generated processed data
released under CC BY 4.0 as described in the repository data licence. The
original workbooks should be deposited with the article or in the immutable
project data archive, rather than duplicated as opaque binary inputs in Git.
