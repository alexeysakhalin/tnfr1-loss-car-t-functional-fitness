# TNFR1 loss and CAR-T functional fitness: analysis repository

This repository contains analysis code and version-controlled data products for
the manuscript “Tumor cell TNFR1 loss attenuates inflammatory responsiveness
and is associated with reduced CAR-T cell functional fitness in an
antigen-retaining in vitro co-culture model.”

The analyses cover tumor-cell bulk RNA-seq, targeted single-cell mRNA profiling,
descriptive DepMap context and exploratory analyses of four published
immune-checkpoint-blockade cohorts. The public cohorts are not CAR-T-treated
cohorts and do not test a TNFR1-dependent clinical mechanism. Transferred C0-C9
scores are bulk expression signatures rather than cell fractions or measured
CAR-T phenotypes in patients.

## Repository structure

| Path | Contents |
|---|---|
| `R/` | Figure-generating R scripts and shared plotting/validation functions |
| `scripts/` | Source acquisition, deterministic table preparation and bulk RNA-seq modelling |
| `validation/` | Independent CheckMate survival recalculation |
| `resources/` | Frozen C0-C9 signatures, the versioned R/05 concordance contract, curated gene sets and identifier mappings |
| `data/experimental/` | Version-controlled experimental analysis tables and provenance |
| `data/analysis/` | Aggregate CheckMate results and the compact DepMap S1B derivative |
| `data/source_manifest.tsv` | Public-source locations, licences, expected filenames, sizes and SHA-256 values |
| `reference_results/` | Aggregate numerical reference results without sample identifiers |
| `docs/` | Analysis details, figure map and complete input inventory |
| `tests/` | Structural and numerical regression tests |

The exact local-only inputs and their required filenames are listed in
[`docs/REPRODUCIBILITY_INPUTS.md`](docs/REPRODUCIBILITY_INPUTS.md).

## Software environment

The automated workflows use Python 3.12 and R 4.4.3. The core and bulk
RNA-seq requirement files pin different dependency versions, so install them
in separate environments from the repository root.

Core validation and open-cohort preparation:

```bash
python -m venv .venv-analysis
source .venv-analysis/bin/activate
python -m pip install --upgrade pip
python -m pip install --requirement requirements.txt
deactivate
```

Bulk RNA-seq reconstruction:

```bash
python -m venv .venv-bulk
source .venv-bulk/bin/activate
python -m pip install --upgrade pip
python -m pip install --requirement requirements-bulk-rnaseq.txt
deactivate
```

`renv.lock` records the R environment used by the automated bulk RNA-seq and
DepMap renderers; it is not a complete lock for the single-cell or clinical
scripts. The targeted single-cell script additionally requires `Matrix`,
`Seurat`, `dplyr`, `tidyr`, `tibble`, `patchwork`, `writexl`,
`data.table`, `ggplot2` and `scales`. Clinical-context scripts declare their
package requirements at the start of each file. Archive `sessionInfo()` with
each final figure run.

## Validation

Run the lightweight repository checks with:

```bash
python -m py_compile scripts/*.py validation/*.py tests/*.py
python -m unittest discover -s tests -v
```

GitHub Actions runs the same contracts, parses all R scripts and independently
rebuilds the bulk RNA-seq and DepMap S1B figures. Figure artifacts include PNG,
600-dpi LZW TIFF, numerical output contracts, source checksums and runtime
provenance.

The targeted single-cell CI also validates both per-cluster marker workbooks
before artifact assembly. The sheet, row-count and OOXML integrity contract is
documented in [`docs/XLSX_OUTPUT_CONTRACT.md`](docs/XLSX_OUTPUT_CONTRACT.md).

## Bulk RNA-seq

The tracked count matrix contains 46,427 Gene_ID features and 24 samples:
WT and TNFR1-KO1, four treatment conditions and three independent experiments
per genotype-condition combination. Source label `T6` maps to manuscript clone
`TNFR1-KO1`. The three records displayed as `TRNAV-CAC` are summed only for the
symbol-level models, yielding 46,425 unique symbols.

Re-estimate the seven primary contrasts and three prespecified
genotype-by-treatment interactions, then render Figures 1B, 1C, 2B and
Supplementary Figure S2D:

```bash
python scripts/run_bulk_rnaseq_pydeseq2.py \
  --input-dir data/experimental/bulk_rnaseq \
  --output-dir results/bulk_rnaseq \
  --include-interaction \
  --n-cpus 2
Rscript R/render_bulk_rnaseq_figures.R
```

PyDESeq2 v0.5.4 fits each primary contrast as a separate six-sample model.
Complete 46,425-symbol result universes, including explicitly non-estimable
rows, are passed to the figure scripts. The within-treatment TNFR1-KO1-versus-WT
contrasts used in Figure 2B are distinct from the formal interaction tests.

The validated complete adapters are also version-controlled under
`data/experimental/bulk_rnaseq/derived/`, so the R figure scripts can be run
directly from a clean clone. The display contracts verify panel order, gene
universe, missingness, thresholds, Venn membership, label selection, image
dimensions and TIFF compression.

## Targeted single-cell mRNA profiling

`R/05_figure_4_AB_suppl_S5A.R` reads the four tracked `CD3+ cells` matrices:
WT (3,743 cells), KO1 (2,475), KO2 (1,831) and the repeated-stimulation dataset
(5,662). Each matrix contains 259 targeted genes. These data must be described
as targeted single-cell mRNA profiling rather than whole-transcriptome
single-cell RNA-seq.

```bash
Rscript R/05_figure_4_AB_suppl_S5A.R
```

Before figure release, review the exported marker profiles against every manual
cluster label and archive the before/after QC counts, marker tables and
`sessionInfo()`.

## Clinical-context analyses

List the required publisher files and prepare the local analysis tables:

```bash
python scripts/fetch_public_sources.py --list
python scripts/fetch_public_sources.py \
  --source imvigor210_processed_package \
  --accept-licensed-public-downloads
Rscript scripts/export_imvigor210_inputs.R \
  --package-tarball data/raw/IMvigor210CoreBiologies_1.0.0.tar.gz \
  --output-dir data/raw
python scripts/prepare_open_cohort_analysis_tables.py --include-checkmate-aggregates
```

The IMvigor210 exporter verifies the official version 1.0.0 package archive,
recreates the clinical and `log2(CPM + 1)` inputs, and publishes them locally
only when their byte sizes and SHA-256 values match `data/source_manifest.tsv`.
Because the package stores a legacy `DESeq` `CountDataSet`, this one export step
requires a compatible legacy R/Bioconductor environment. Use
`--verify-only --output-dir data/raw` to check existing exports without loading
the package.

Then run the clinical-context scripts in manuscript order:

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

The CheckMate analyses use the 181 nivolumab-treated RNA-profiled tumors
(CM-009, 16; CM-010, 45; CM-025, 120), not the pooled 311-sample dataset.
`OS_CNSR` and `PFS_CNSR` are retained as supplied: 1 denotes an observed event
and 0 denotes censoring. The C6 overall-survival comparison is nominal and
exploratory; its multiplicity-adjusted result is not significant.

The independent recalculation accepts the official combined Braun workbook or
the checksum-verified split S1/S4 files:

```bash
python validation/recalculate_checkmate_survival.py \
  --braun-workbook data/raw/41591_2020_839_MOESM2_ESM.xlsx \
  --output-dir results/checkmate_validation
```

Patient-level clinical and expression data remain local and are excluded by
`.gitignore`. This preserves the source access and redistribution conditions;
only aggregate, non-identifying validation results are version-controlled.
The IMvigor210 package archive and its two sample-level CSV exports are not
included in the GitHub release or the DOI-backed project-data archive.

The project source workbooks used to construct the tracked experimental tables
are deposited separately in Zenodo under the reserved version DOI
[`10.5281/zenodo.19707614`](https://doi.org/10.5281/zenodo.19707614). The record
contains project-generated source workbooks, not third-party clinical cohorts,
DepMap source files or the IMvigor210 package and exports. The DOI will resolve
publicly after the record is published.

## DepMap Supplementary Figure S1B

The tracked compact derivative contains 1,591 eligible cell-line models joined
by `ModelID`. At the prespecified threshold of <0.5 log2(TPM+1), RIPK3 is below
threshold in 1,003/1,591 (63.0%), NLRP3 in 1,172/1,591 (73.7%) and both in
749/1,591 (47.1%). These are descriptive frequencies in the checksum-pinned
source pair, not cancer prevalence or evidence of pathway competence.

Render the panel from a clean clone with:

```bash
Rscript R/11_supplementary_1B.R
```

The expression matrix and `Model.csv` are a confirmed DepMap Public 25Q2
source pair downloaded from the portal's **All Data** page. Exact checksums and
the machine-readable release lock are documented in
[`docs/DEPMAP_S1B.md`](docs/DEPMAP_S1B.md). This portal-hosted release has no
release-specific Figshare DOI.

## Figure map and experimental scope

[`docs/FIGURE_MAP.md`](docs/FIGURE_MAP.md) maps each manuscript panel to its
script, inputs and output filename. Supplementary Figure S3 is a flow-cytometry
experiment and is outside this transcriptomic repository; its final analysis
requires the raw FCS files, compensation/gating record, denominator, replicate
map and prespecified contrasts.

## Licence and citation

Repository code is available under the MIT License. Project-derived processed
data and figures are available under CC BY 4.0; third-party data retain their
source terms. See `LICENSE`, `DATA_LICENSE.md` and `CITATION.cff`.

The clinical-context analyses use processed data accompanying Mariathasan et
al. (2018), Liu et al. (2019), Braun et al. (2020) and Ravi et al. (2023).
Exact citations, source URLs, accessions, licences, expected sizes and SHA-256
checksums are recorded in `data/source_manifest.tsv`.
