# Preparation and verification scripts

| Entry point | Purpose |
|---|---|
| `prepare_experimental_analysis_tables.py` | Build tracked experimental tables from author source workbooks |
| `prepare_bulk_rnaseq_counts.py` | Prepare the bulk count matrix and sample design |
| `run_bulk_rnaseq_pydeseq2.py` | Fit primary and interaction bulk RNA-seq models |
| `prepare_open_cohort_analysis_tables.py` | Prepare local published-cohort inputs and aggregate CheckMate results |
| `prepare_depmap_s1b.py` | Build the compact DepMap S1B derivative |
| `fetch_public_sources.py` | List and retrieve explicitly permitted source files |
| `export_imvigor210_inputs.R` | Recreate version-locked local IMvigor210 package exports |
| `verify_*` | Enforce numerical, workbook, provenance and release contracts |

Canonical filenames, checksums, licences and redistribution policies are in
[`data/source_manifest.tsv`](../data/source_manifest.tsv). Publisher cohort
files and patient-level derivatives remain local.
