# Marker-workbook output contract

`R/05_figure_4_AB_suppl_S5A.R` writes two publication-supporting workbooks:

- `Supplementary_Table_S5_top_markers_per_cluster.xlsx` contains the reviewed
  historical C0-C9 marker ranking; and
- `Supplementary_Table_TCR_top_markers_per_cluster.xlsx` contains the
  independently clustered repeated-stimulation C0-C5 marker ranking.

Each workbook has `All_clusters_top20` first, followed by one `Cluster_<n>`
worksheet in numeric cluster order. The summary has 20 rows per cluster and
each cluster worksheet has exactly 20 data rows. Column order is inherited
unchanged from the pinned Seurat marker table.

The workbooks are written directly from these tables with `writexl 2.0.0`.
The targeted single-cell workflow then runs
`scripts/verify_xlsx_workbook.py` before copying any workbook into the uploaded
artifact. The verifier checks ZIP integrity, resolves every internal OOXML
relationship, rejects references to absent parts, requires worksheet dimensions
to match populated cells, and enforces the expected sheet order and table
extents. A malformed workbook therefore fails CI rather than being repaired or
published. The release gate additionally opens each workbook in normal mode
with the repository-pinned `openpyxl 3.1.5` reader and verifies the same sheet
sequence and extents.

The descriptive Figure 4A-B panels retain all QC-passing C0-C10 cells. C10 is
reported neutrally as a small cytokine/IFN-response-high cluster and remains
outside the historical C0-C9 frozen/transferred signature mapping; it is not
classified as a contaminant. The S5 marker workbook remains C0-C9 because its
purpose is to validate that historical mapping. The complete C0-C10 marker
table is retained separately as
`Supplementary_Table_S5_all_cluster_markers_including_C10.tsv.gz`.
