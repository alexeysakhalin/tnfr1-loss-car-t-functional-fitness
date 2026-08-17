# Automated validation

The workflows in this directory are release checks for the complete analysis
repository:

| Workflow | Scope |
|---|---|
| `validate.yml` | Python contracts, source listing and R syntax |
| `bulk-rnaseq.yml` | Bulk RNA-seq model reconstruction and manuscript panels |
| `cohort-inputs-targeted-singlecell.yml` | Licensed cohort-input verification and targeted single-cell outputs |
| `depmap-s1b.yml` | DepMap Public 25Q2 Supplementary Figure S1B |

Raw publisher cohorts and patient-level derivatives are never uploaded as
artifacts. Release artifacts contain figures, aggregate results, checksums and
runtime provenance only.

The targeted single-cell workflow uses an explicit filename allowlist. Its C6
projection outputs are limited to the rendered exploratory figure, aggregate
cluster summaries and frozen-signature gene coverage. Cell barcodes, cell-level
scores, expression matrices and serialized Seurat objects are not staged.
The versioned annotation manifest records manual post-clustering marker
interpretation and explicitly states that no reference-atlas classifier or
cross-dataset label transfer was used.
