# Experimental analysis tables

This directory contains compact, non-identifying, author-generated tables used
by the experimental figure scripts. The eight source Excel workbooks are not
used directly by R and are not committed here. Their canonical names, byte
sizes and SHA-256 values are recorded in `experimental_data_manifest.tsv`.
The manifest uses stable canonical deposition names. Timestamp suffixes added
by file transfer are removed from those names; the byte-level source identity
is preserved by SHA-256. The longer supplied IFN filename containing
`no_zero_pvals` is normalized because two p-values are zero from numerical
underflow. Deposit the unchanged source workbooks under the manifest's stable
canonical names; SHA-256, rather than a transient upload suffix, defines their
identity.

To rebuild these tables from the source workbooks:

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

The significant-gene table supports the thresholded overlap in Figure 1C and
checks of reported gene-level effects. It does not support a conventional
volcano plot or a complete tested-gene enrichment universe. Figure 1B remains
locked until full, unfiltered DESeq2 outputs are provided under the schema
enforced by `R/03_figure_1_B_C.R`.

These tables are processed DE results, not raw-count inputs. The repository
can reproduce the plotted thresholds from them, but re-estimation of the
underlying DE models additionally requires the count matrix, design/contrast
code and locked DESeq2/PyDESeq2 environment in the article-data archive.

## Matched T6-versus-WT RNA-seq

The four `hela_t6_vs_wt_*_differential_expression.tsv.gz` files contain the
complete condition-specific tables from the matched workbook. The associated
design contains 24 samples: HeLa and T6, four conditions, three biological
replicates per cell-line/condition combination. These files support Figure 2B
and Supplementary Figure S2D.

`T6` is the internal cell-line identifier encoded in the workbook. The file
does not establish whether T6 corresponds to the manuscript label KO1 or KO2.
The code therefore preserves `T6`; a signed sample key or author confirmation
is required before replacing that identifier in the manuscript or figure.

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
