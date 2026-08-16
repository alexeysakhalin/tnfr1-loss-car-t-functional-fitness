# Analysis data and input policy

Version-controlled experimental tables and non-identifying aggregate results
are stored under `data/experimental/` and `data/analysis/`. Publisher cohort
matrices and sample-level clinical tables are not distributed. They are
downloaded or supplied locally, verified against
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
| `data/source_manifest.tsv` | Source URLs, citations, accessions, licences, canonical filenames/checksums and verification notes | yes |
| `data/experimental/` | Non-identifying project-derived transcriptomic analysis tables, checksums and QC | yes |
| `data/depmap/raw/` | Checksum-pinned DepMap expression ZIP and `Model.csv` used for S1B | no |
| `data/analysis/depmap_s1b_eligible_models.tsv.gz` | Eligible DepMap S1B models, expression, flags and quadrants | yes |
| `data/analysis/depmap_s1b_preparation_qc.json` | Source-selection, join and denominator QC for S1B | yes |
| `data/analysis/depmap_s1b_source_provenance.json` | Raw-source hashes and release-pair status for S1B | yes |
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
- `41591_2019_654_MOESM4_ESM.xlsx` and either
  `41591_2019_654_MOESM3_ESM.txt` or the checksum-pinned
  `Liu2019_NatureMedicine_metastatic_melanoma_antiPD1_expression_matrix.csv`
  (Liu clinical table and TPM matrix);
- `IMvigor210_clinical.csv` and `IMvigor210_expression_log2CPM.csv`, recreated
  locally from the checksum-pinned `IMvigor210CoreBiologies` 1.0.0 package as
  described below; the clinical file has an exact byte gate, while expression
  uses the versioned semantic contract;
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
downloader. The preparation script consumes two locally generated project
exports: an exact-checksum clinical data frame and a semantically verified
log2-CPM expression matrix. Recreate them from the official package archive
with:

```bash
python scripts/fetch_public_sources.py \
  --source imvigor210_processed_package \
  --accept-licensed-public-downloads
Rscript scripts/export_imvigor210_inputs.R \
  --package-tarball data/raw/IMvigor210CoreBiologies_1.0.0.tar.gz \
  --output-dir data/raw
python scripts/verify_imvigor210_expression.py \
  --input data/raw/IMvigor210_expression_log2CPM.csv
```

The exporter first verifies the 122,127,298-byte package archive against
SHA-256 `cfdd3176d7b34de5b04fb9416bfd2b20fa4b6e238aaad5f20b048a34329ea178`.
It loads the package `cds` object, writes `Biobase::pData(cds)` as the clinical
table, and calculates `log2(edgeR::cpm(DESeq::counts(cds)) + 1)` on the complete
31,286-feature count matrix. Every cell must agree with the independent direct
formula `log2(sweep(counts, 2, colSums(counts), "/") * 1e6 + 1)` and with the
CSV write/readback within `5e-13 + 5e-14 * abs(expected)`; ordered identifiers,
dimensions, finite/nonnegative values and the structural-zero mask must also
match. It requires the pinned legacy R 4.0 / Bioconductor 3.11 environment.

The package archive and clinical CSV retain exact size/SHA-256 gates. The
expression CSV instead must match the fixed6 digest in
`resources/IMvigor210_expression_semantic_contract_v1.json`; fixed7 and fixed8
hashes are diagnostics. Its manifest size and SHA-256 identify the canonical
historical rendering but are not acceptance gates. The streaming verifier uses
only Python's standard library and writes no identifiers or expression values
to its report. Existing local files can be checked with:

```bash
Rscript scripts/export_imvigor210_inputs.R \
  --verify-only --output-dir data/raw
```

`scripts/prepare_open_cohort_analysis_tables.py` repeats the expression
semantic verification before reading the matrix and converts every value to
the required fixed6 `ROUND_HALF_UP` representation before identifier mapping,
aggregation, within-sample ranking or output. Therefore sub-six-decimal text
differences cannot alter a rank or another downstream result.

The official package archive and both patient/sample-level CSV exports remain
local-only. They are excluded from Git and from the DOI-backed code/data
archive, including Zenodo; readers obtain the licensed package from its
recorded source and run the exporter. Entrez identifiers are subsequently
mapped with the checksum-pinned HGNC snapshot. Any separately redistributed
derivative must retain attribution, the CC BY 3.0 licence and a change notice.

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
any sample-level derivative is committed. The primary published-cohort analysis is
restricted to the 181 nivolumab-treated RNA samples (CM-009: 16; CM-010: 45;
CM-025: 120). Only aggregate coefficients, confidence intervals, multiplicity
adjustments and QC counts are version-controlled.


## DepMap Supplementary Figure S1B

The DepMap source is not patient-cohort data and is not combined with the four
clinical cohorts. The complete expression archive and source `Model.csv`
remain local-only. Exact filenames, byte sizes, SHA-256 values, the metadata
MD5 and redistribution policy are recorded in `source_manifest.tsv`.

`scripts/prepare_depmap_s1b.py` selects 1,684
`is_default_entry=True` expression records, one per unique `ModelID`, joins
them to `Model.csv`, and retains 1,591 records annotated as `Cell Line` with a
non-empty OncoTree primary disease other than `Non-Cancerous`. `TissueOrigin`
is empty in the supplied metadata and is not used as a filter. The compact
tracked derivative contains only the model/profile identifiers, disease label,
two expression values, threshold flags and quadrant assignment.

`R/11_supplementary_1B.R` reads that tracked derivative, QC, provenance and
statistics contract, so a clean clone can render the panel without the full
raw files. The authors confirmed that both source files were downloaded from
the DepMap Portal **All Data** page for DepMap Public 25Q2. Their checksums,
confirmed release-pair status and `same_release_pair=true` assertion are
machine-locked. The complete contract, fixed counts and wording restrictions
are documented in `docs/DEPMAP_S1B.md`.
