# TNFR1 loss and CAR-T functional fitness: analysis repository

This repository contains the analysis code accompanying the manuscript
“Tumor cell TNFR1 loss attenuates inflammatory responsiveness and reduces
CAR-T cell functional fitness in an antigen-retaining in vitro co-culture
model.”

The code follows the three evidence levels used in the revised manuscript:

1. the tumor-cell phenotype, including cytokine-induced apoptotic signalling
   and ICAM1 induction after TNFR1 loss;
2. the later in-vitro CAR-T phenotype after serial exposure to
   TNFR1-deficient targets; and
3. exploratory clinical-context analyses of four published
   immune-checkpoint-blockade cohorts.

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
| `data/experimental/` | Non-identifying experimental analysis tables and QC manifest |
| `tests/` | Repository-contract and numerical-regression tests |

Sample-level clinical and expression data are not version-controlled. Source
files and locally generated sample-level tables are excluded by `.gitignore`.
This avoids redistributing the Braun CheckMate supplement, which does not carry
an open redistribution licence, and keeps one consistent policy across all
four cohorts.

Source Excel workbooks should be kept unchanged in the article data archive. R
scripts consume compact canonical TSV.gz files rather than opaque binary
workbooks. Patient-level public-cohort inputs remain local-only.

## Software environment

Python 3.12 and R 4.4 are used by the automated checks; the bulk RNA-seq
release workflow fixes Python 3.12.13. Create a Python
environment and install the pinned validation dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install --requirement requirements.txt
```

The final article release must additionally include the `renv.lock` and
`sessionInfo()` produced by the clean R/Seurat run. Package versions are not
guessed before that run is completed.

The bulk RNA-seq model is isolated in a separately pinned Python environment:

```bash
python -m pip install --requirement requirements-bulk-rnaseq.txt
```

## Experimental workflow

The canonical experimental tables are already version-controlled. Their
source-workbook hashes and transformation QC are documented in
[`data/experimental/README.md`](data/experimental/README.md). Authors can
rebuild them from the eight source workbooks with:

```bash
python scripts/prepare_experimental_analysis_tables.py --input-dir /path/to/workbooks
```

The original bulk RNA-seq workbook has been converted into a deterministic
Gene_ID-level integer-count matrix with 46,427 features and 24 samples. The
sample manifest records WT and TNFR1-KO1 cells, four treatment conditions and
three independent experiments per combination. The internal source clone name
`T6` is mapped to manuscript clone `TNFR1-KO1`; the only repeated display
symbol, `TRNAV-CAC`, is summed over its three source Gene_ID records only when
the symbol-level models are fitted. Source and canonical-table checksums are
recorded under [`data/experimental/bulk_rnaseq/`](data/experimental/bulk_rnaseq/).

Re-estimate all seven primary contrasts and the three prespecified
genotype-by-treatment interaction contrasts, then generate the experimental
bulk-transcriptomic panels in manuscript order:

```bash
python scripts/run_bulk_rnaseq_pydeseq2.py \
  --input-dir data/experimental/bulk_rnaseq \
  --output-dir results/bulk_rnaseq \
  --include-interaction \
  --n-cpus 2
Rscript R/03_figure_1_B_C.R
Rscript R/04_figure_2_B_suppl_S2D.R
Rscript R/05_figure_4_AB_suppl_S5A.R
```

PyDESeq2 v0.5.4 fits each primary contrast as a separate six-sample model:
three WT cytokine-versus-untreated comparisons and four TNFR1-KO1-versus-WT
comparisons within treatment strata. Complete 46,425-symbol result universes,
including non-significant and explicitly non-estimable rows, are supplied to
`R/03` and `R/04`. Positive log2 fold changes follow the first-named
numerator-versus-reference direction, and the exported Wald statistic is
checked against `log2FoldChange/lfcSE`. Optional full-factorial interaction
outputs are kept separate from the within-stratum contrasts used in Figure 2B.
The earlier significant-only Excel exports are retained only as provenance and
are not used to draw the final volcano plots.

The two complete figure-input adapters and the three complete interaction
tables from the validated release run are version-controlled under
`data/experimental/bulk_rnaseq/derived/`. Consequently, R03 and R04 can render
the bulk panels directly from a clean clone. The count-level workflow rebuilds
these files and verifies exact structural, missingness, threshold-category,
Venn-membership and prespecified interaction-gene outcomes before rendering
from the committed release adapters. Raw cross-platform numeric deltas remain
available in the comparison report; the displayed volcano-label genes have an
additional absolute log2-fold-change reproducibility bound of `0.001`.

The `Rebuild bulk RNA-seq figures` GitHub Actions workflow performs the pinned
count-to-DE run, checks all complete result exports, renders 600-dpi PNG/TIFF
panels with R 4.4.3 and stores both result and figure artifacts together with
the comparison report and runtime provenance.

The targeted single-cell inputs use the `CD3+ cells` sheets: WT, 3,743 cells;
KO1, 2,475; KO2, 1,831; and the repeated-stimulation dataset, 5,662. These are
259-gene targeted count matrices, not whole-transcriptome single-cell RNA-seq.
The clean release run must archive before/after QC counts, nonempty marker
tables and author-reviewed marker-to-label mappings; equality of numeric
cluster identifiers alone is not accepted as annotation validation.

## Clinical-context workflow

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

## Release-locked analyses

Supplementary Figure S1B is intentionally release-locked. The supplied
expression archive has the profile-level `is_default_entry` schema associated
with DepMap Public 25Q2, but the archive itself does not identify its release.
The analysis therefore still requires the same-release `Model.csv`, the
official release/source URL, download date, checksums and one-default-profile-
per-human-cell-line filtering. A release DOI is recorded when available but is
not mandatory because current DepMap releases are distributed directly from
the portal rather than Figshare. The repository does not retain the earlier
unverified DepMap percentages.

Supplementary Figure S3 is not generated by the transcriptomic workflow. It
requires raw FCS files, the compensation/gating record, a defined live-cell
denominator, replicate and donor mapping, and the exact statistical contrasts
before the panel and legend can be treated as final.

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

The lightweight workflow checks syntax and data contracts. A release is made
only after a clean R run regenerates the mapped panels, the figures are
compared with the source tables, and the resulting `sessionInfo()` is archived.

## Licence

Repository code is available under the MIT License. Author-generated processed
data and figures are available under CC BY 4.0; third-party data retain their
source terms. See `LICENSE` and `DATA_LICENSE.md`.

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
