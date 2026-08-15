# Supplementary Figure S1B: DepMap source and analysis contract

Supplementary Figure S1B is a descriptive analysis of two checksum-pinned
DepMap files: a profile-level expression matrix and `Model.csv`. The code
selects one default RNA-expression profile per `ModelID`, joins annotations,
and applies the population rule defined below.

## Publication status

The numerical contract is complete and the expression matrix is fixed as
DepMap Public 25Q2. The release identity of the separately checksum-pinned
`Model.csv` must still be confirmed from the authors' original DepMap download
record before a same-release claim is made. Its checksum is independently
reproduced by a frozen copy labelled DepMap Public 25Q3, but an official
historical DepMap inventory exposing the 25Q2 `Model.csv` checksum was not
available for verification. This does not prove that the 25Q2 and 25Q3
metadata files differ; it means only that the supplied metadata release is
currently unverified.

There are two publication-safe resolutions:

1. Obtain `Model.csv` from the same verified DepMap release as the expression
   matrix, rerun the checksum/count contract, and update the fixed values if
   that file differs.
2. Use the currently supplied pair and identify the expression and metadata
   releases separately in Methods, the figure legend and provenance. A
   cross-release metadata join is acceptable as a descriptive annotation step
   only when it is disclosed explicitly.

Do not silently label both files as DepMap Public 25Q2. DepMap Public 25Q2 was
distributed through the DepMap portal without a release-specific Figshare DOI;
a DOI must not be invented.

## Checksum-pinned inputs

| Input | Repository-local filename | Size (bytes) | Digest |
|---|---|---:|---|
| Supplied expression ZIP | `data/depmap/raw/OmicsExpressionTPMLogp1HumanProteinCodingGenes.zip` | 249,426,032 | SHA-256 `c44524c48e20f8c5c1263eb23cd55df77ceda62cfb5246babbe22cecc90c3da0` |
| CSV member | `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` | 538,420,733 | SHA-256 `90bfdbe5c44cbb8f822e655ba7f179f3033933116285b6b2f85153b2d3d17c75` |
| Supplied metadata | `data/depmap/raw/Model.csv` | 699,474 | SHA-256 `9dbb9de8805696c1345816ab07edd23fb4fd95e117739f3c5c3b1cf062c1233b`; MD5 `af4472ab734ea3aec974d992b504c7e5` |

The full files remain local-only. The expression ZIP exceeds GitHub's normal
per-file limit, and the full annotation table is unnecessary for figure
reproduction. `scripts/prepare_depmap_s1b.py` validates both raw sources and
creates four deterministic, version-controlled artifacts:

| Tracked artifact | Purpose | SHA-256 |
|---|---|---|
| `data/analysis/depmap_s1b_eligible_models.tsv.gz` | 1,591 eligible rows, expression values, threshold flags and quadrant | `368ad92b085a722d3984a5355bea3109d8e5a2b29ffe563c3fd284cf8970f354` |
| `data/analysis/depmap_s1b_preparation_qc.json` | Source, join and denominator QC | `8d58a7113ca08c7ffd8297e26f3a3a29693e61e119186365c1fb04531d988b79` |
| `data/analysis/depmap_s1b_source_provenance.json` | Source digests, population rule and release-status lock | `e21b63f90835e80e71c388b13e167b214773cd6a988aa43a3aaac8cd65745242` |
| `reference_results/depmap_s1b_statistics.csv` | Six fixed counts and percentages | `60585de1dd22e879220d7d9da89d6cd762685b1b11af1d60b64627490e80e990` |

The compressed TSV is sorted by `ModelID` and written with gzip modification
time zero and no original filename. Repeated preparation from unchanged source
files therefore produces identical bytes. The source-provenance artifact fixes
the expression release as 25Q2, the metadata release identity as unverified,
`release_pair_status` as `unverified`, and `same_release_pair` as null. Release
labels are not accepted as render-time input; changing this state requires an
intentional, reviewed repository update after source confirmation.

## Population rule and audit counts

The expression matrix has 1,739 profile rows: 1,684 records with
`is_default_entry == True` and 55 non-default records. The 1,684 default rows
contain 1,684 unique `ModelID` and 1,684 unique `ProfileID` values.

The supplied `Model.csv` has 2,132 unique `ModelID` rows. It contains 2,108
`Cell Line` and 24 `Organoid` records; 151 rows have
`OncotreePrimaryDisease == "Non-Cancerous"`. `TissueOrigin` is present but
empty in all 2,132 rows. It is therefore recorded as unavailable and is not a
required field or a filter.

All 1,684 default-expression `ModelID` values match one metadata row. All are
annotated `ModelType == "Cell Line"`; 93 are annotated
`OncotreePrimaryDisease == "Non-Cancerous"`. The analysis retains records for
which:

- `is_default_entry == True`;
- `ModelType == "Cell Line"`;
- `OncotreePrimaryDisease` is non-empty and is not `"Non-Cancerous"`; and
- both RIPK3 and NLRP3 values are finite.

This yields 1,591 eligible DepMap cell-line models. The denominator must not
be called "human" because no species filter was possible from the supplied
metadata. The most exact manuscript label is:

> DepMap cell-line models with a non-missing OncoTree primary-disease label
> other than `Non-Cancerous`.

## Fixed descriptive results

The prespecified threshold is `<0.5 log2(TPM+1)`.

| Metric | n / 1,591 | Percent |
|---|---:|---:|
| RIPK3 below threshold | 1,003 | 63.0% |
| NLRP3 below threshold | 1,172 | 73.7% |
| Both below threshold | 749 | 47.1% |
| RIPK3 below threshold only | 254 | 16.0% |
| NLRP3 below threshold only | 423 | 26.6% |
| Neither below threshold | 165 | 10.4% |

These are frequencies in the checksum-pinned source pair, not prevalence
estimates for cancers or patients. Use "below the prespecified threshold"
rather than "absent" or an unqualified "low". Transcript abundance does not
establish necroptotic or pyroptotic pathway competence.

## Reproduction

A clean clone renders the manuscript panel directly from tracked data:

```bash
Rscript R/11_supplementary_1B.R
```

This writes one manuscript-ready scatter plot as PNG and 600-dpi LZW TIFF,
plus runtime provenance and `sessionInfo()`, under
`results/supplementary_S1B/`. The renderer checks all four tracked artifact
sizes, the row schema, stored flags, quadrant assignments, denominator and six
reported counts before plotting. Repository tests independently pin and verify
the SHA-256 value of every tracked artifact.

To regenerate the tracked derivative from the full sources, place the
unchanged files at the local-only paths above and run:

```bash
python scripts/prepare_depmap_s1b.py
```

The preparation script validates both raw-file sizes and digests, the sole ZIP
member and its digest, all source schemas and counts, one-default-profile-per-
model selection, metadata coverage, population filters and threshold results.
It then rewrites the TSV, QC, source-provenance and statistics artifacts
deterministically. Release labels are not accepted as command-line overrides:
until confirmed, they cannot be injected into the locked provenance by a
rendering command.

## Manuscript wording

Methods should state the release of each input separately unless same-release
provenance is verified. A concise template is:

> We selected records with `is_default_entry=True` from the checksum-pinned
> DepMap Public 25Q2 profile-level log2(TPM+1) matrix and joined them by
> `ModelID` to a separately checksum-pinned `Model.csv` [insert its verified
> release, or state that the metadata release could not be verified]. We retained
> records annotated as `ModelType="Cell Line"` with a non-missing
> `OncotreePrimaryDisease` other than `Non-Cancerous`. `TissueOrigin` was empty
> in the supplied metadata and was not used. RIPK3 and NLRP3 were classified as
> below threshold at <0.5 log2(TPM+1).

Results and the legend may state:

> Among 1,591 eligible DepMap cell-line models, RIPK3 was below the
> prespecified threshold in 1,003 (63.0%), NLRP3 in 1,172 (73.7%), and both in
> 749 (47.1%). These descriptive frequencies do not establish pathway
> competence or prevalence in human cancers.

In the Discussion, use this panel only as broad descriptive context for the
experimental tumor-cell and CAR-T phenotypes. It should not be presented as
mechanistic or clinical validation.

## Citation

Cite each verified DepMap release in the form requested by the portal, for
example `DepMap, Broad (2025). DepMap Public 25Q3. Dataset. depmap.org`, and
cite the portal at <https://depmap.org/portal/>. Also cite:

Arafeh R, Shibue T, Dempster JM, et al. The present and future of the Cancer
Dependency Map. *Nature Reviews Cancer*. 2025;25:59-73.
<https://doi.org/10.1038/s41568-024-00763-x>
