# Figure and analysis scripts

The numbered R scripts follow manuscript figure order where practical.

| Scripts | Scope |
|---|---|
| `00_*` | Shared data loading and validation helpers |
| `01`, `02`, `06`-`10`, `12` | Exploratory published-cohort panels in Figure 5 and Supplementary Figure S6 |
| `03`, `04` | Experimental bulk RNA-seq panels in Figures 1-2 and Supplementary Figure S2 |
| `05` | Targeted single-cell panels in Figures 4-5 and Supplementary Figure S5 |
| `11` | DepMap Supplementary Figure S1B |
| `render_bulk_rnaseq_figures.R` | Release renderer for bulk RNA-seq panels |

Exact inputs and output filenames for every reproduced panel are listed in
[`docs/FIGURE_MAP.md`](../docs/FIGURE_MAP.md). These scripts do not reproduce
wet-laboratory flow-cytometry panels.
