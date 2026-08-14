# TNFR1 loss and CAR-T functional fitness: analysis repository

This repository contains the analysis code accompanying the manuscript
“Tumor cell TNFR1 loss attenuates inflammatory responsiveness and reduces
CAR-T cell functional fitness in an antigen-retaining in vitro co-culture
model.”

The code is organized around two distinct evidence layers:

1. experimental bulk and targeted single-cell analyses from the in-vitro
   HeLa–CAR-T system; and
2. exploratory clinical-context analyses of four published immune-checkpoint
   blockade cohorts.

The public cohorts are not CAR-T-treated cohorts and are not used to validate a
TNFR1-dependent clinical mechanism. Transferred C0–C9 scores are bulk
expression signatures, not cell fractions or measured CAR-T phenotypes in
patients.

## Repository contents

| Path | Purpose |
|---|---|
| `R/` | Figure-generating R scripts and shared validation helpers |
| `scripts/` | Checksum-verified source acquisition and deterministic cohort preparation |
| `validation/` | Independent Python recalculation of the primary CheckMate survival analyses |
| `resources/` | Frozen C0–C9 marker definitions, curated display genes and identifier mappings |
| `reference_results/` | Aggregate validation results without sample identifiers |
| `data/source_manifest.tsv` | Source URL, DOI/accession, licence, expected filename, size and SHA-256 |
| `tests/` | Repository-contract and numerical-regression tests |

Sample-level clinical and expression data are not version-controlled. Source
files and locally generated sample-level tables are excluded by `.gitignore`.
This avoids redistributing the Braun CheckMate supplement, which does not carry
an open redistribution licence, and keeps one consistent policy across all
four cohorts.

## Clinical-context workflow

Python 3.12 and R 4.4 are used by the automated checks. Create a Python
environment and install the pinned validation dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install --requirement requirements.txt
```

Review the source inventory and obtain the required publisher files:

```bash
python scripts/fetch_public_sources.py --list
```

The detailed acquisition commands, licences, local filenames and schemas are
documented in [`data/README.md`](data/README.md). Each input is rejected unless
its byte size and SHA-256 match `data/source_manifest.tsv`.

After the verified inputs have been placed under `data/raw/`, prepare the local
analysis tables and aggregate CheckMate model results:

```bash
python scripts/prepare_open_cohort_analysis_tables.py --include-checkmate-aggregates
```

No SQLite database or monolithic R data object is required. The clinical
figures can then be generated in manuscript order:

```bash
Rscript R/01_validate_analysis_tables.R
Rscript R/02_supplementary_S6A_S6C.R
Rscript R/06_figure_5_C.R
Rscript R/07_figure_5D.R
Rscript R/08_figure_5E_S6D_S6E.R
Rscript R/09_figure_5F_S6G.R
Rscript R/10_figure_5G.R
Rscript R/12_supplementary_S6B_S6F.R
```

The independent CheckMate recalculation accepts either the official combined
Braun workbook or the verified split S1/S4 files:

```bash
python validation/recalculate_checkmate_survival.py \
  --braun-workbook data/raw/41591_2020_839_MOESM2_ESM.xlsx \
  --output-dir results/checkmate_validation
```

## Prespecified CheckMate population and event coding

The publisher workbook contains 311 RNA-profiled tumors: 181 from the
nivolumab arm and 130 from the everolimus arm. Clinical-context analyses in the
manuscript use only the 181 nivolumab-treated tumors (CM-009, 16; CM-010, 45;
CM-025, 120). `OS_CNSR` and `PFS_CNSR` are retained as supplied: 1 denotes an
observed event and 0 denotes censoring.

The cross-cohort bulk T-cell score uses the five genes measured in every
cohort: `CD2`, `CD3D`, `CD3E`, `CD8A` and `CD8B`.

The C6 median-split overall-survival comparison is nominal and exploratory:
trial-stratified Cox HR 1.504 (95% CI 1.050–2.154), Cox *p*=0.0259 and
BH-adjusted Cox *p*=0.166 across ten signatures. The descriptive log-rank
*p*=0.0246 has BH-adjusted *p*=0.208. No C6 progression-free-survival
association was detected. These results do not support a predictive biomarker
claim.

## Experimental inputs

Scripts for Figures 1, 2, 4 and the repeated-stimulation component of Figure 5
require project RNA-seq and targeted single-cell input files that are not part
of the patient-cohort archive. Their required local paths are listed in
[`docs/FIGURE_MAP.md`](docs/FIGURE_MAP.md). The scripts fail explicitly when a
required input or package is missing.

Supplementary Figure S1B is intentionally release-locked: it must not be used
until the exact DepMap Public release, Figshare DOI, download date, checksums
and one-default-profile-per-human-cell-line filtering are supplied. The
repository does not retain the earlier unverified DepMap percentages.

## Validation

Run the lightweight checks locally with:

```bash
python -m py_compile scripts/*.py validation/*.py tests/*.py
python -m unittest discover -s tests -v
```

GitHub Actions independently parses every R script and runs the Python
contract tests on each pull request. Aggregate reference results are retained
to detect changes in cohort composition, event coding, the frozen 20-gene
membership of each C0-C9 signature or multiplicity adjustment.

## Source studies

The clinical-context analyses use processed data accompanying:

- Mariathasan et al., *Nature* (2018),
  [doi:10.1038/nature25501](https://doi.org/10.1038/nature25501);
- Liu et al., *Nature Medicine* (2019),
  [doi:10.1038/s41591-019-0654-5](https://doi.org/10.1038/s41591-019-0654-5);
- Braun et al., *Nature Medicine* (2020),
  [doi:10.1038/s41591-020-0839-y](https://doi.org/10.1038/s41591-020-0839-y);
- Ravi et al., *Nature Genetics* (2023),
  [doi:10.1038/s41588-023-01355-5](https://doi.org/10.1038/s41588-023-01355-5).

Exact publisher URLs, accessions, licences, expected byte sizes and SHA-256
checksums are recorded in `data/source_manifest.tsv`.
