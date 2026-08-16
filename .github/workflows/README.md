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
