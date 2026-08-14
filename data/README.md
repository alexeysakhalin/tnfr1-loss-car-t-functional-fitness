# Public-cohort data workflow

This repository does not distribute source matrices or sample-level clinical
tables. They are downloaded or supplied locally, verified against
`source_manifest.tsv`, and transformed by deterministic scripts. Only source
metadata, transformation code, gene-set definitions, aggregate quality-control
results and aggregate model coefficients are version-controlled.

This policy applies to all four clinical cohorts. It prevents an open software
repository from becoming an uncontrolled copy of patient-level data and avoids
mixing publisher supplements, controlled-access raw data and project-derived
exports in one directory.

## Local directory layout

| Path | Contents | Version-controlled |
|---|---|---|
| `data/source_manifest.tsv` | Source URLs, citations, accessions, licences, expected filenames and SHA-256 values | yes |
| `data/experimental/` | Non-identifying project-derived transcriptomic analysis tables, checksums and QC | yes |
| `resources/analysis_gene_identifiers.tsv` | Frozen identifiers for the 181 analysis genes, derived from the declared HGNC snapshot | yes |
| `resources/hgnc_20260814_gene_identifiers.tsv.gz` | Frozen full HGNC mapping used before transcriptome-wide ranking | yes |
| `data/raw/` | Publisher files and package exports | no |
| `data/analysis/open_cohort_selected_expression.tsv.gz` | Combined sample-level selected expression for IMvigor210, SU2C-MARK and Liu | no |
| `data/analysis/open_cohort_sample_metadata.tsv.gz` | Combined minimal metadata for the three open cohorts | no |
| `data/processed/checkmate_selected_expression.tsv.gz` | Local selected expression for all 311 CheckMate RNA samples | no |
| `data/processed/checkmate_sample_metadata.tsv.gz` | Local CheckMate metadata for all 311 RNA samples | no |
| `data/analysis/checkmate_c6_global_gene_models.tsv.gz` | Aggregate gene-level coefficients from the nivolumab CheckMate analysis | yes |
| `data/analysis/checkmate_c6_group_balance.tsv` | Aggregate residualized-group balance statistics | yes |

The experimental tables have a separate provenance and schema guide in
`data/experimental/README.md`. They are not patient-cohort data and are
version-controlled under the repository data licence.

## Reproduction

From the repository root:

```bash
python scripts/fetch_public_sources.py --list
python scripts/fetch_public_sources.py --source braun_checkmate_supplement
```

The Braun workbook is placed under `data/raw/` only after its size and SHA-256
match the manifest. The identifier maps required by the analysis are frozen in
`resources/`; the manifest records the full-source HGNC snapshot. Other public processed
sources can be downloaded from the pinned URLs in the manifest. Controlled raw
sequencing accessions are provenance records, not downloader inputs.

Place these verified files in `data/raw/` before running the preparation step:

- `41588_2023_1355_MOESM3_ESM.xlsx` (SU2C-MARK Supplementary Tables 1 and 13);
- `41591_2019_654_MOESM4_ESM.xlsx` and
  `41591_2019_654_MOESM3_ESM.txt` (Liu clinical table and TPM matrix);
- `IMvigor210_clinical.csv` and `IMvigor210_expression_log2CPM.csv`, the
  checksum-pinned project exports from `IMvigor210CoreBiologies` 1.0.0;
- `41591_2020_839_MOESM2_ESM.xlsx` (downloaded by the script for local
  CheckMate calculations).

Then run once:

```bash
python scripts/prepare_open_cohort_analysis_tables.py --include-checkmate-aggregates
```

This writes the two combined `open_cohort_*` tables and a local QC record, then
writes the two local CheckMate sample-level tables and recalculates the
aggregate CheckMate C6 gene models from the official Braun workbook. Omit
`--include-checkmate-aggregates` when only the three open-cohort tables are
needed. No sample-level result is staged for Git.

## Table schemas

Both selected-expression tables contain:

| Column | Definition |
|---|---|
| `cohort_id` | Stable cohort identifier used by the analysis scripts |
| `sample_id` | Source pseudonymous RNA/sample identifier |
| `gene_symbol` | Requested analysis symbol after deterministic identifier mapping |
| `expr_value` | Source expression value; no cross-cohort rescaling is applied here |
| `expression_unit` | Cohort-specific unit (`log2CPM`, `TPM`, or `normalized_expression`) |
| `rank_percentile` | Within-sample percentile after all mapped source features are aggregated to unique gene symbols; calculated before analysis-gene selection |

Both metadata tables contain `cohort_id`, `sample_id`, `trial_id`,
`treatment_arm`, OS/PFS time and event columns, and `tumor_purity`. Fields not
used or not harmonized for an open cohort remain empty rather than being
silently inferred. CheckMate event coding is retained exactly as documented by
the source: 1 is an observed event and 0 is censored. These tables remain local
even when a publisher licence permits redistribution.

The cross-cohort bulk T-cell score uses the five genes present in every source
matrix: `CD2`, `CD3D`, `CD3E`, `CD8A` and `CD8B`. `TRAC` is not available in
IMvigor210, Liu or CheckMate and is therefore not silently treated as zero.

The aggregate CheckMate model table contains one row per tested gene and no
sample identifiers. Its coefficient compares residualized C6-high with
residualized C6-low tumors in the 181 nivolumab-treated RNA samples, adjusting
for the bulk T-cell score and trial. Genes used to calculate either score are
excluded from the model family; `BH_p` is adjusted over all remaining fitted
genes.

## Cohort-specific decisions

### IMvigor210

The public processed package is `IMvigor210CoreBiologies` 1.0.0 (CC BY 3.0).
Raw sequencing data at EGA are controlled access and are not used by the
downloader. The preparation script consumes two checksum-pinned project
exports: the clinical data frame and the log2-CPM expression matrix. The exact
historical commands that produced these CSV files from the package were not
present in the supplied repository; reproduction therefore begins from the
two verified exports. They should be deposited in the versioned data archive
before publication. Entrez identifiers are mapped with the checksum-pinned
HGNC snapshot. Any redistributed derivatives must retain attribution and state
that they were reformatted and gene-filtered.

### SU2C-MARK

The public Supplementary Table workbook is used directly. The clinical header
is row 3 of `Table_S1_Clinical_Annotations`; the RNA TPM header is row 3 of
`Table_S13_RNA_TPM`. Reading the workbook directly avoids locale conversion of
decimal points and 27 spreadsheet-formatted gene labels found in an earlier
CSV export. Stable Ensembl identifiers are stripped of version suffixes and
mapped through the frozen HGNC table when the displayed symbol is unusable.
All 152 RNA samples are retained; the 136 `Keep` and 16 `Flag` source labels
are counted explicitly in the local preparation QC record.

### Liu melanoma

The public table and TPM supplements are used directly. Only clinical rows
whose first field matches `^Patient[0-9]+$` are data records; footer notes are
never interpreted as patients. All 121 TPM samples must match the clinical
table. The source is attributed under CC BY 4.0 and the local derivative is
described as reformatted and gene-filtered.

### Braun CheckMate

The official combined Supplementary Table workbook is downloaded and verified
locally. Checksum-verified split S1/S4 copies are accepted as equivalent local
inputs. The source has no open redistribution licence. Neither the workbook nor
any sample-level derivative is committed. The primary clinical-context analysis is
restricted to the 181 nivolumab-treated RNA samples (CM-009: 16; CM-010: 45;
CM-025: 120). Only aggregate coefficients, confidence intervals, multiplicity
adjustments and QC counts are version-controlled.

## Why SQLite is not used

The analysis consumes four fixed cohort matrices and a small number of explicit
transformations. Two explicit compressed TSV contracts plus a source manifest
make every transformation inspectable without a database schema or hidden
import state. Full matrices remain in their publisher formats; selected local
tables are regenerated when gene sets or source versions change.
