# Experimental bulk RNA-seq

This directory contains the complete integer read-count matrix used for the
HeLa CD19 bulk RNA-seq analyses. The source workbook was supplied by Macrogen
and is recorded by byte size and SHA-256 in `source_manifest.tsv`.

## Files

| File | Contents |
|---|---|
| `gene_counts.tsv.gz` | One unique source `Gene_ID` per row and 24 integer-count sample columns |
| `gene_annotations.tsv.gz` | Source annotation corresponding one-to-one with `gene_counts.tsv.gz` |
| `gene_symbol_membership.tsv.gz` | Deterministic mapping from display symbols to source `Gene_ID` values |
| `sample_metadata.tsv` | Stable sample identifiers, genotype, treatment and replicate |
| `source_manifest.tsv` | Provenance and checksum of the source workbook |
| `source_qc.json` | Structural, count-integrity, duplicate-symbol and per-sample QC |
| `SHA256SUMS` | Checksums of the canonical tables and QC record |
| `derived/` | Complete validated figure adapters, interaction tables and release-run metadata |

The source contains 46,427 unique `Gene_ID` rows and 46,425 unique gene
symbols. `TRNAV-CAC` is the only repeated symbol: it represents Gene_ID
107985614, 107985615 and 107985753. Canonical counts remain at `Gene_ID`
resolution. The manuscript analysis is performed at symbol level by summing
counts over Gene_ID rows before model fitting; the rule is deterministic and
changes no other symbol.

`T6` is the source clone name for the manuscript TNFR1-KO1 clone. `TI` denotes
combined TNF and IFNγ treatment. These aliases are retained in
`sample_metadata.tsv` and are not inferred from filenames during analysis.
Treatments were 50 ng/mL TNF, 50 ng/mL IFNγ or both for 48 hours. Replicates
are recorded as independent experiments. Whether WT and KO1 replicate numbers
represent paired experimental batches was not confirmed, so no pairing term is
used or implied.

## Rebuild the canonical tables

With the checksum-matched source workbook available locally:

```bash
python scripts/prepare_bulk_rnaseq_counts.py \
  --input Expression_Profile.GRCh38.gene.xlsx \
  --output-dir data/experimental/bulk_rnaseq
```

The preparation step rejects a changed checksum, altered sheet structure,
missing values, negative or fractional counts, duplicate Gene_ID values and
any duplicate symbol other than the documented `TRNAV-CAC` mapping. The local
workbook basename may differ: source identity is established by SHA-256, while
the exact filename originally received and the canonical archive filename are
frozen in `source_manifest.tsv`. That provenance row is author-maintained and
is not overwritten by the table-conversion command.

## Differential-expression analysis

Install the pinned optional environment and run:

```bash
python -m pip install --requirement requirements-bulk-rnaseq.txt
python scripts/run_bulk_rnaseq_pydeseq2.py \
  --input-dir data/experimental/bulk_rnaseq \
  --output-dir results/bulk_rnaseq \
  --include-interaction \
  --n-cpus 2
```

This generates three WT cytokine-versus-control contrasts and four
TNFR1-KO1-versus-WT contrasts. Every comparison uses a separate six-sample
two-group fit, matching the normalization scope of the legacy tables. All
46,425 symbol-level features are written to every output. No filtering by
base mean, fold change, p-value or adjusted p-value is applied to the exported
tables; expected missing statistics are retained and labelled.

The `--include-interaction` option fits the prespecified
`genotype + treatment + genotype:treatment` model and export three formal
difference-in-differences contrasts. The interaction results are distinct from
the within-stratum contrasts used for the manuscript figures.

The committed `derived/` snapshot was produced with Python 3.12.13,
PyDESeq2 0.5.4 and the exact environment in
`derived/environment.freeze.txt`. It contains complete, unfiltered
46,425-symbol universes rather than selected significant rows.
`derived/SHA256SUMS` freezes this release snapshot; the original count matrix and model
code remain the source of truth.
