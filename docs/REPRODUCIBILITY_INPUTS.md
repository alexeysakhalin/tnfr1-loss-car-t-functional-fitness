# Reproducibility input inventory

This inventory distinguishes files that are included in the repository from
source files that must be supplied locally. Exact hashes and sizes are enforced
by `data/experimental/experimental_data_manifest.tsv`,
`data/experimental/bulk_rnaseq/source_manifest.tsv` and
`data/source_manifest.tsv`.

## Analysis-ready inputs included in the repository

| Analysis | Repository input |
|---|---|
| Bulk RNA-seq count model | `data/experimental/bulk_rnaseq/gene_counts.tsv.gz`; `gene_annotations.tsv.gz`; `gene_symbol_membership.tsv.gz`; `sample_metadata.tsv` |
| Figure 1B/1C direct render | `data/experimental/bulk_rnaseq/derived/figure_1b_1c_wt_cytokine_contrasts.unfiltered.tsv.gz` |
| Figure 2B/S2D direct render | `data/experimental/bulk_rnaseq/derived/figure_2b_s2d_tnfr1_ko1_vs_wt_matched_treatments.unfiltered.tsv.gz` |
| Interaction analyses | three `data/experimental/bulk_rnaseq/derived/interaction_*.unfiltered.tsv.gz` files |
| Targeted single-cell panels | `data/experimental/singlecell/WT_targeted_counts.tsv.gz`; `KO1_targeted_counts.tsv.gz`; `KO2_targeted_counts.tsv.gz`; `TCR_targeted_counts.tsv.gz` |
| DepMap S1B direct render | `data/analysis/depmap_s1b_eligible_models.tsv.gz`; preparation QC; source provenance; `reference_results/depmap_s1b_statistics.csv` |
| CheckMate Figure 5F/S6G | `data/analysis/checkmate_c6_global_gene_models.tsv.gz`; `data/analysis/checkmate_c6_group_balance.tsv` |
| Fixed analysis resources | `resources/CAR_T_state_signatures.csv`; `resources/CAR_T_state_signature_concordance_v1.csv`; `resources/Figure_5F_curated_gene_sets.csv`; identifier maps |

## Project source workbooks kept outside GitHub

Extract the Zenodo source-workbook archive and pass its root directory to the
preparation script. The script accepts either a flat directory or the deposited
`bulk_rnaseq/` and `targeted_single_cell/` subdirectories, searches recursively,
and requires exactly one match for every filename pattern. Filenames must match
exactly.

| Required filename | Source content | Used by |
|---|---|---|
| `Expression_Profile.GRCh38.gene.xlsx` | Macrogen GRCh38 gene-level integer counts; 46,427 features x 24 samples | `scripts/prepare_bulk_rnaseq_counts.py` |
| `Differential_Expression_TNFa_vs_control_final_filtered.xlsx` | historical TNF-versus-control DE export | provenance/QC cross-check in `scripts/prepare_experimental_analysis_tables.py` |
| `Differential_Expression_IFN_vs_control_filtered.xlsx` | historical IFN-gamma-versus-control DE export | provenance/QC cross-check in `scripts/prepare_experimental_analysis_tables.py` |
| `Differential_Expression_TI_vs_control_final_filtered.xlsx` | historical TNF+IFN-gamma-versus-control DE export | provenance/QC cross-check in `scripts/prepare_experimental_analysis_tables.py` |
| `PyDESeq2_T6_vs_Hela_matched_treatments.xlsx` | historical matched T6-versus-WT four-condition DE export | provenance/QC cross-check; `T6` = `TNFR1-KO1` |
| `WT_new.xlsx` | WT CD3-positive targeted count matrix | `scripts/prepare_experimental_analysis_tables.py` |
| `KO1_new.xlsx` | TNFR1-KO1 CD3-positive targeted count matrix | same script |
| `KO2_new.xlsx` | TNFR1-KO2 CD3-positive targeted count matrix | same script |
| `TCR.xlsx` | repeated-stimulation targeted count workbook | same script; `CD3+ cells` sheet, 5,662 cells |

The historical DE workbooks are not the gene universe for the final volcano
plots. Figures 1B/1C and 2B/S2D use the complete count-level PyDESeq2 results.

Rebuild canonical experimental tables with:

```bash
python scripts/prepare_experimental_analysis_tables.py --input-dir /path/to/workbooks
python scripts/prepare_bulk_rnaseq_counts.py \
  --input /path/to/workbooks/bulk_rnaseq/Expression_Profile.GRCh38.gene.xlsx \
  --output-dir data/experimental/bulk_rnaseq
```

For a flat local directory, omit `bulk_rnaseq/` from the second command.

## Public-cohort inputs kept outside GitHub

Place verified publisher files under `data/raw/` using the exact filenames
below. `scripts/fetch_public_sources.py --list` reports source URLs, licences,
canonical sizes and SHA-256 values. These are exact gates except where the
IMvigor210 expression row explicitly identifies its byte hash as provenance
and points to the required semantic contract.

| Required filename | Source |
|---|---|
| `41591_2020_839_MOESM2_ESM.xlsx` | Braun et al. CheckMate ccRCC supplement; automatic checksum-verified download supported |
| `41588_2023_1355_MOESM3_ESM.xlsx` | Ravi et al. SU2C-MARK supplement |
| `41591_2019_654_MOESM4_ESM.xlsx` | Liu et al. clinical supplement |
| one of `41591_2019_654_MOESM3_ESM.txt` or `Liu2019_NatureMedicine_metastatic_melanoma_antiPD1_expression_matrix.csv` | Liu et al. expression supplement, either as the official text file or the checksum-pinned CSV export accepted by the preparer |
| `IMvigor210CoreBiologies_1.0.0.tar.gz` | official CC BY 3.0 processed package; supplied to `scripts/export_imvigor210_inputs.R`, not to the cohort-preparation script |
| `IMvigor210_clinical.csv` | exact-checksum clinical export generated locally by `scripts/export_imvigor210_inputs.R` |
| `IMvigor210_expression_log2CPM.csv` | semantically pinned `log2(CPM + 1)` export generated locally by the same script; fixed6 contract is required, canonical byte hash is provenance only |

These files contain publisher or sample-level data and are not redistributed
through the repository. The preparation script writes local-only harmonized
tables under `data/analysis/` and `data/processed/`.

Download the official IMvigor210 package archive from the source URL in
`data/source_manifest.tsv`, retain its canonical filename, and run:

```bash
python scripts/fetch_public_sources.py \
  --source imvigor210_processed_package \
  --accept-licensed-public-downloads
Rscript scripts/export_imvigor210_inputs.R \
  --package-tarball data/raw/IMvigor210CoreBiologies_1.0.0.tar.gz \
  --output-dir data/raw
python scripts/verify_imvigor210_expression.py \
  --input data/raw/IMvigor210_expression_log2CPM.csv
Rscript scripts/export_imvigor210_inputs.R \
  --verify-only --output-dir data/raw
```

The downloader verifies the package archive against the manifest. The exporter
verifies it again before reading it. It retains exact byte gates for that archive
and the clinical CSV. For expression it checks all cells against the independent
direct CPM formula, repeats the all-cell check after CSV write/readback, and
requires the fixed6 semantic digest defined in
`resources/IMvigor210_expression_semantic_contract_v1.json`. Fixed7/fixed8 and
the canonical file SHA-256 are provenance diagnostics only. The streaming
Python verifier retains no expression matrix and its report contains dimensions
and hashes, never sample/feature identifiers or expression values. The cohort
preparer then canonicalizes every accepted cell to fixed6 `ROUND_HALF_UP`
before mapping and ranking; fixed6-compatible CSV variants consequently produce
the same analysis matrix rather than merely passing the same input check.
The historical-input impact audit is committed as
`resources/IMvigor210_fixed6_canonicalization_impact_v1.json`: all 61,944
selected-expression keys and every rank percentile were unchanged; 59,954
printed expression values changed by at most `4.999993601373376e-7`. Figure 5C
scores, S6 score orderings, Spearman results and significance/selection
thresholds were unchanged. Density-rendering bytes may change on regeneration,
but no reported interpretation changes.

Loading the legacy `CountDataSet` requires the pinned R 4.0 / Bioconductor 3.11
environment with `DESeq`, `Biobase` and `edgeR`; current analyses use R 4.4.3.
The final R command checks already generated files without loading the package
and invokes the semantic verifier. The package archive and both CSV files remain
local-only and are not included in the GitHub release or DOI-backed archive,
including Zenodo.

## DepMap source files kept outside GitHub

Place the following files under `data/depmap/raw/` only when regenerating the
tracked compact derivative:

| Required filename | Source identity |
|---|---|
| `OmicsExpressionTPMLogp1HumanProteinCodingGenes.zip` | DepMap Public 25Q2 expression archive |
| `Model.csv` | checksum-pinned DepMap Public 25Q2 model metadata downloaded from the same **All Data** release page as the expression archive |

Run:

```bash
python scripts/prepare_depmap_s1b.py
Rscript R/11_supplementary_1B.R
```

The compact derivative is already included, so the full DepMap files are not
required to render Supplementary Figure S1B from a clean clone.

## Separate experimental source for Supplementary Figure S3

Supplementary Figure S3 is not produced by the transcriptomic scripts. Its
reproducibility package consists of raw FCS files, the compensation and gating
workspace, the live-cell denominator, E:T ratio and timing, replicate/donor
mapping and the prespecified statistical contrasts.
