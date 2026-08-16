# Repository contract tests

The test suite guards the released analysis against changes in:

- input structure, checksums and declared licences;
- bulk RNA-seq contrasts, thresholds and manuscript result universes;
- targeted single-cell marker/concordance and workbook contracts;
- published-cohort and DepMap aggregate results;
- workflow privacy rules and artifact contents.

Run from the repository root:

```bash
python -m unittest discover -s tests -v
```

The tests are release safeguards; they do not replace biological interpretation
or review of wet-laboratory replicate structure.
