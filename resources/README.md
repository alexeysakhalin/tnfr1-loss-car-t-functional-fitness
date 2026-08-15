# Frozen analysis resources

`CAR_T_state_signatures.csv` contains ten fixed C0-C9 expression signatures,
with 20 positively ranked genes per computational cluster. The membership was
frozen before the revised public-cohort scoring and is checksum-locked by the
repository tests.

These are top-ranked genes from an exploratory cell-level marker analysis,
not ten sets of biological-replicate differential-expression findings. No
adjusted-p-value filter is applied when the fixed 20-gene membership is read.
In the frozen table, 198 of 200 rows have adjusted p-values below 0.05; KLRC1
in C0 and CHI3L2 in C9 do not. This does not change membership, but it prevents
the signatures from being described as 200 statistically significant markers.

The single-cell experiment has one targeted count matrix per conditioning
label. Marker p-values and cluster fractions are therefore descriptive at the
cell/computational level and do not provide biological-replicate inference or
absolute viable-cell phenotype frequencies.

`Figure_5F_curated_gene_sets.csv` defines the prespecified display categories
for the exploratory aggregate CheckMate gene-model plot. Identifier maps are
checksum-pinned snapshots used during cohort preparation.
