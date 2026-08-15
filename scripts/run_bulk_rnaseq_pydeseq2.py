#!/usr/bin/env python3
"""Run the prespecified bulk RNA-seq differential-expression analyses."""

from __future__ import annotations

import argparse
import atexit
import gzip
import hashlib
import io
import json
import platform
import shutil
import tempfile
from importlib.metadata import distributions, version
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


REQUIRED_PYDESEQ2_VERSION = "0.5.4"
TREATMENTS = ("control", "IFNG", "TNF", "TNF_IFNG")
RESULT_COLUMNS = (
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
)
WT_FIGURE_COLUMNS = (
    "condition",
    "gene_symbol",
    "base_mean",
    "log2_fold_change_treatment_vs_untreated",
    "lfc_se",
    "wald_statistic_treatment_vs_untreated",
    "p_value",
    "adjusted_p_value",
)
KO_FIGURE_COLUMNS = (
    "condition",
    "gene_symbol",
    "base_mean",
    "log2_fold_change_ko1_vs_wt",
    "lfc_se",
    "wald_statistic_ko1_vs_wt",
    "p_value",
    "adjusted_p_value",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_dataframe_gzip(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as text_handle:
                frame.to_csv(
                    text_handle,
                    sep="\t",
                    index=False,
                    lineterminator="\n",
                    na_rep="NA",
                    float_format="%.17g",
                )


def load_inputs(input_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    counts_path = input_dir / "gene_counts.tsv.gz"
    annotations_path = input_dir / "gene_annotations.tsv.gz"
    metadata_path = input_dir / "sample_metadata.tsv"
    for path in (counts_path, annotations_path, metadata_path):
        if not path.exists():
            raise FileNotFoundError(path)

    counts = pd.read_csv(counts_path, sep="\t", dtype={"gene_id": "string"})
    annotations = pd.read_csv(
        annotations_path,
        sep="\t",
        dtype="string",
        keep_default_na=False,
    )
    metadata = pd.read_csv(metadata_path, sep="\t", dtype="string")

    if counts["gene_id"].duplicated().any():
        raise ValueError("gene_counts.tsv.gz contains duplicate gene_id values")
    if annotations["gene_id"].duplicated().any():
        raise ValueError("gene_annotations.tsv.gz contains duplicate gene_id values")
    if set(counts["gene_id"]) != set(annotations["gene_id"]):
        raise ValueError("Count and annotation Gene_ID universes differ")
    if metadata["sample_id"].duplicated().any():
        raise ValueError("sample_metadata.tsv contains duplicate sample_id values")

    sample_ids = metadata["sample_id"].tolist()
    if counts.columns.tolist() != ["gene_id", *sample_ids]:
        raise ValueError("Count columns must exactly match sample_metadata.tsv order")
    numeric_counts = counts[sample_ids].apply(pd.to_numeric, errors="raise")
    if not np.isfinite(numeric_counts.to_numpy()).all():
        raise ValueError("Count matrix contains non-finite values")
    if (numeric_counts.to_numpy() < 0).any():
        raise ValueError("Count matrix contains negative values")
    if not np.equal(numeric_counts.to_numpy(), np.floor(numeric_counts.to_numpy())).all():
        raise ValueError("Count matrix contains non-integer values")
    counts[sample_ids] = numeric_counts.astype(np.int64)

    metadata["replicate"] = pd.to_numeric(metadata["replicate"], errors="raise").astype(int)
    if set(metadata["genotype"]) != {"WT", "TNFR1_KO1"}:
        raise ValueError("Expected WT and TNFR1_KO1 genotypes")
    if set(metadata["treatment"]) != set(TREATMENTS):
        raise ValueError(f"Expected treatments {TREATMENTS!r}")
    group_sizes = metadata.groupby(["genotype", "treatment"], observed=True).size()
    if len(group_sizes) != 8 or not (group_sizes == 3).all():
        raise ValueError("Expected three replicates in every genotype-by-treatment group")
    return counts, annotations, metadata


def aggregate_to_gene_symbols(
    counts: pd.DataFrame,
    annotations: pd.DataFrame,
    sample_ids: Sequence[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    joined = counts.merge(
        annotations[["gene_id", "gene_symbol"]],
        on="gene_id",
        how="left",
        validate="one_to_one",
    )
    if joined["gene_symbol"].isna().any() or (joined["gene_symbol"] == "").any():
        raise ValueError("Every Gene_ID must have a non-empty gene_symbol")

    symbol_counts = (
        joined.groupby("gene_symbol", sort=True, observed=True)[list(sample_ids)]
        .sum()
        .astype(np.int64)
    )
    membership = (
        joined.groupby("gene_symbol", sort=True, observed=True)["gene_id"]
        .agg(
            n_source_gene_ids="size",
            source_gene_ids=lambda values: ";".join(
                sorted(values.astype(str), key=lambda value: int(value))
            ),
        )
    )
    if symbol_counts.shape != (46_425, 24):
        raise ValueError(
            f"Expected a 46,425 by 24 symbol count matrix; observed {symbol_counts.shape}"
        )
    duplicate_membership = membership[membership["n_source_gene_ids"] > 1]
    if duplicate_membership.index.tolist() != ["TRNAV-CAC"]:
        raise ValueError("Unexpected repeated gene symbol after canonical aggregation")
    if duplicate_membership.loc["TRNAV-CAC", "source_gene_ids"] != (
        "107985614;107985615;107985753"
    ):
        raise ValueError("TRNAV-CAC Gene_ID membership differs from the source contract")
    return symbol_counts, membership


def make_metadata(metadata: pd.DataFrame, sample_ids: Sequence[str]) -> pd.DataFrame:
    selected = metadata.set_index("sample_id").loc[list(sample_ids), ["genotype", "treatment"]]
    selected = selected.copy()
    genotype_categories = [
        value for value in ("WT", "TNFR1_KO1") if value in set(selected["genotype"])
    ]
    treatment_categories = [
        value for value in TREATMENTS if value in set(selected["treatment"])
    ]
    selected["genotype"] = pd.Categorical(
        selected["genotype"], categories=genotype_categories, ordered=True
    )
    selected["treatment"] = pd.Categorical(
        selected["treatment"], categories=treatment_categories, ordered=True
    )
    return selected


def fit_model(
    symbol_counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design: str,
    n_cpus: int,
) -> tuple[DeseqDataSet, pd.DataFrame]:
    sample_ids = metadata.index.tolist()
    sample_by_gene = symbol_counts[sample_ids].T
    nonzero = sample_by_gene.sum(axis=0) > 0
    model_counts = sample_by_gene.loc[:, nonzero]
    dds = DeseqDataSet(
        counts=model_counts,
        metadata=metadata,
        design=design,
        refit_cooks=True,
        n_cpus=n_cpus,
        quiet=True,
        low_memory=True,
    )
    dds.deseq2()
    return dds, sample_by_gene


def complete_results(
    stats: DeseqStats,
    sample_by_gene: pd.DataFrame,
    membership: pd.DataFrame,
) -> pd.DataFrame:
    stats.summary()
    observed = stats.results_df.copy()
    observed.index = observed.index.astype(str)
    missing_columns = set(RESULT_COLUMNS).difference(observed.columns)
    if missing_columns:
        raise ValueError(f"PyDESeq2 result columns are missing: {sorted(missing_columns)!r}")

    universe = pd.DataFrame(index=sample_by_gene.columns.astype(str))
    universe.index.name = "Gene_Symbol"
    complete = universe.join(observed[list(RESULT_COLUMNS)], how="left")
    complete = complete.join(membership, how="left")
    total_count = sample_by_gene.sum(axis=0).reindex(complete.index).astype(np.int64)
    all_zero = total_count.eq(0)
    # PyDESeq2 omits all-zero features from the fitted object. Their arithmetic
    # mean is nevertheless defined and equals zero; only inferential fields are
    # non-estimable. Restoring baseMean=0 also preserves standard DESeq2 output
    # semantics for downstream universe and plotting checks.
    complete.loc[all_zero, "baseMean"] = 0.0
    if complete.loc[all_zero, [
        "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"
    ]].notna().any().any():
        raise ValueError("All-zero features unexpectedly contain test statistics")
    complete.insert(0, "total_count_in_model_subset", total_count)
    complete.insert(
        0,
        "analysis_status",
        np.where(
            total_count.eq(0),
            "all_zero_in_subset",
            np.where(complete["baseMean"].notna(), "modelled", "not_returned_by_engine"),
        ),
    )
    complete = complete.reset_index()
    complete = complete[
        [
            "Gene_Symbol",
            "source_gene_ids",
            "n_source_gene_ids",
            "total_count_in_model_subset",
            "analysis_status",
            *RESULT_COLUMNS,
        ]
    ]
    if len(complete) != 46_425 or complete["Gene_Symbol"].duplicated().any():
        raise ValueError("Every unfiltered output must contain 46,425 unique symbols")
    zero_rows = complete["analysis_status"].eq("all_zero_in_subset")
    if not complete.loc[zero_rows, "baseMean"].eq(0.0).all():
        raise ValueError("All-zero features must be exported with baseMean equal to zero")

    finite = complete[["log2FoldChange", "lfcSE", "stat"]].notna().all(axis=1)
    finite &= complete["lfcSE"].ne(0)
    if finite.any():
        expected_stat = complete.loc[finite, "log2FoldChange"] / complete.loc[finite, "lfcSE"]
        if not np.allclose(
            complete.loc[finite, "stat"], expected_stat, rtol=1e-6, atol=1e-8
        ):
            raise ValueError("Wald statistic direction is inconsistent with log2FoldChange")
    for column in ("pvalue", "padj"):
        values = complete[column].dropna()
        if not values.between(0, 1, inclusive="both").all():
            raise ValueError(f"{column} contains values outside [0, 1]")
    return complete.sort_values("Gene_Symbol", kind="stable").reset_index(drop=True)


def run_contrast(
    dds: DeseqDataSet,
    sample_by_gene: pd.DataFrame,
    membership: pd.DataFrame,
    contrast: list[str] | np.ndarray,
    output_path: Path,
    n_cpus: int,
) -> tuple[dict[str, Any], pd.DataFrame]:
    stats = DeseqStats(
        dds,
        contrast=contrast,
        alpha=0.05,
        cooks_filter=True,
        independent_filter=True,
        n_cpus=n_cpus,
        quiet=True,
    )
    complete = complete_results(stats, sample_by_gene, membership)
    write_dataframe_gzip(complete, output_path)
    status_counts = complete["analysis_status"].value_counts().sort_index().to_dict()
    return {
        "output_file": output_path.name,
        "rows": len(complete),
        "analysis_status_counts": status_counts,
        "sha256": sha256_file(output_path),
        "contrast": contrast.tolist() if isinstance(contrast, np.ndarray) else contrast,
    }, complete


def figure_adapter(
    complete: pd.DataFrame,
    treatment: str,
    comparison: str,
) -> pd.DataFrame:
    adapter = complete[["Gene_Symbol", *RESULT_COLUMNS]].copy()
    adapter.insert(0, "condition", treatment)
    common_mapping = {
        "Gene_Symbol": "gene_symbol",
        "baseMean": "base_mean",
        "lfcSE": "lfc_se",
        "pvalue": "p_value",
        "padj": "adjusted_p_value",
    }
    if comparison == "treatment_vs_untreated":
        mapping = {
            **common_mapping,
            "log2FoldChange": "log2_fold_change_treatment_vs_untreated",
            "stat": "wald_statistic_treatment_vs_untreated",
        }
        return adapter.rename(columns=mapping)[list(WT_FIGURE_COLUMNS)]
    if comparison == "ko1_vs_wt":
        mapping = {
            **common_mapping,
            "log2FoldChange": "log2_fold_change_ko1_vs_wt",
            "stat": "wald_statistic_ko1_vs_wt",
        }
        return adapter.rename(columns=mapping)[list(KO_FIGURE_COLUMNS)]
    raise ValueError(f"Unknown figure comparison: {comparison}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("data/experimental/bulk_rnaseq"),
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results/bulk_rnaseq")
    )
    parser.add_argument("--n-cpus", type=int, default=1)
    parser.add_argument(
        "--include-interaction",
        action="store_true",
        help="Fit genotype-by-treatment difference-in-differences contrasts",
    )
    args = parser.parse_args()

    observed_version = version("pydeseq2")
    if observed_version != REQUIRED_PYDESEQ2_VERSION:
        raise RuntimeError(
            f"This workflow requires pydeseq2 {REQUIRED_PYDESEQ2_VERSION}; "
            f"observed {observed_version}"
        )
    if args.n_cpus < 1:
        raise ValueError("--n-cpus must be at least 1")
    requested_output_dir = args.output_dir
    requested_output_dir.parent.mkdir(parents=True, exist_ok=True)
    if requested_output_dir.exists():
        if not requested_output_dir.is_dir():
            raise FileExistsError(
                f"Output path exists and is not a directory: {requested_output_dir}"
            )
        if any(requested_output_dir.iterdir()):
            raise FileExistsError(
                "Output directory is not empty; move or remove the previous run first: "
                f"{requested_output_dir}"
            )
        requested_output_dir.rmdir()
    staging_output_dir = Path(
        tempfile.mkdtemp(
            prefix=f".{requested_output_dir.name}.",
            suffix=".tmp",
            dir=requested_output_dir.parent,
        )
    )

    def cleanup_staging() -> None:
        shutil.rmtree(staging_output_dir, ignore_errors=True)

    atexit.register(cleanup_staging)
    args.output_dir = staging_output_dir

    counts, annotations, metadata = load_inputs(args.input_dir)
    sample_ids = metadata["sample_id"].tolist()
    symbol_counts, membership = aggregate_to_gene_symbols(counts, annotations, sample_ids)
    manifest_rows: list[dict[str, Any]] = []
    wt_adapters: list[pd.DataFrame] = []
    ko_adapters: list[pd.DataFrame] = []
    figure_input_dir = args.output_dir / "figure_inputs"

    # Separate six-sample fits retain the pairwise normalization and dispersion
    # strategy used for the manuscript WT cytokine contrasts.
    for treatment in TREATMENTS[1:]:
        selected_ids = metadata.loc[
            metadata["genotype"].eq("WT")
            & metadata["treatment"].isin(["control", treatment]),
            "sample_id",
        ].tolist()
        model_metadata = make_metadata(metadata, selected_ids)
        dds, sample_by_gene = fit_model(
            symbol_counts, model_metadata, "~treatment", args.n_cpus
        )
        output_path = args.output_dir / f"WT_{treatment}_vs_control.unfiltered.tsv.gz"
        record, complete = run_contrast(
            dds,
            sample_by_gene,
            membership,
            ["treatment", treatment, "control"],
            output_path,
            args.n_cpus,
        )
        adapter = figure_adapter(complete, treatment, "treatment_vs_untreated")
        wt_adapters.append(adapter)
        write_dataframe_gzip(
            adapter,
            figure_input_dir
            / f"figure_1b_1c_wt_{treatment.lower()}_vs_control.unfiltered.tsv.gz",
        )
        manifest_rows.append(
            {
                "analysis_id": f"WT_{treatment}_vs_control",
                "model": "~treatment",
                "sample_subset": f"WT: control and {treatment}",
                "numerator": f"WT {treatment}",
                "denominator": "WT control",
                "positive_log2fc": f"higher in WT {treatment}",
                **record,
            }
        )

    # These separate six-sample fits reproduce the legacy within-stratum
    # TNFR1-KO1-versus-WT tables used for the manuscript figures.
    for treatment in TREATMENTS:
        selected_ids = metadata.loc[metadata["treatment"].eq(treatment), "sample_id"].tolist()
        model_metadata = make_metadata(metadata, selected_ids)
        dds, sample_by_gene = fit_model(
            symbol_counts, model_metadata, "~genotype", args.n_cpus
        )
        output_path = (
            args.output_dir
            / f"TNFR1_KO1_vs_WT_{treatment}.unfiltered.tsv.gz"
        )
        record, complete = run_contrast(
            dds,
            sample_by_gene,
            membership,
            ["genotype", "TNFR1_KO1", "WT"],
            output_path,
            args.n_cpus,
        )
        adapter = figure_adapter(complete, treatment, "ko1_vs_wt")
        ko_adapters.append(adapter)
        write_dataframe_gzip(
            adapter,
            figure_input_dir
            / (
                "figure_2b_s2d_tnfr1_ko1_vs_wt_"
                f"{treatment.lower()}.unfiltered.tsv.gz"
            ),
        )
        manifest_rows.append(
            {
                "analysis_id": f"TNFR1_KO1_vs_WT_{treatment}",
                "model": "~genotype",
                "sample_subset": f"{treatment}: WT and TNFR1_KO1",
                "numerator": f"TNFR1_KO1 {treatment}",
                "denominator": f"WT {treatment}",
                "positive_log2fc": f"higher in TNFR1_KO1 under {treatment}",
                **record,
            }
        )

    if args.include_interaction:
        full_metadata = make_metadata(metadata, sample_ids)
        dds, sample_by_gene = fit_model(
            symbol_counts,
            full_metadata,
            "~genotype + treatment + genotype:treatment",
            args.n_cpus,
        )
        wt_control = dds.cond(genotype="WT", treatment="control")
        ko_control = dds.cond(genotype="TNFR1_KO1", treatment="control")
        for treatment in TREATMENTS[1:]:
            wt_treated = dds.cond(genotype="WT", treatment=treatment)
            ko_treated = dds.cond(genotype="TNFR1_KO1", treatment=treatment)
            # PyDESeq2 0.5.4 delegates condition-vector construction to
            # formulaic-contrasts, which returns a pandas Series in the locked
            # environment. DeseqStats requires a numeric ndarray for a custom
            # contrast; a Series would otherwise be interpreted as a
            # three-element categorical contrast.
            contrast = np.asarray(
                ko_treated - wt_treated - ko_control + wt_control,
                dtype=float,
            )
            output_path = (
                args.output_dir
                / (
                    "interaction_TNFR1_KO1_vs_WT_"
                    f"{treatment}_vs_control.unfiltered.tsv.gz"
                )
            )
            record, _ = run_contrast(
                dds,
                sample_by_gene,
                membership,
                contrast,
                output_path,
                args.n_cpus,
            )
            manifest_rows.append(
                {
                    "analysis_id": (
                        f"interaction_TNFR1_KO1_vs_WT_{treatment}_vs_control"
                    ),
                    "model": "~genotype + treatment + genotype:treatment",
                    "sample_subset": "all 24 samples",
                    "numerator": (
                        f"(TNFR1_KO1-WT) under {treatment} minus "
                        "(TNFR1_KO1-WT) under control"
                    ),
                    "denominator": "difference-in-differences null of zero",
                    "positive_log2fc": (
                        f"larger KO1-vs-WT difference after {treatment} than at control"
                    ),
                    **record,
                }
            )

    wt_combined = pd.concat(wt_adapters, ignore_index=True)
    write_dataframe_gzip(
        wt_combined,
        figure_input_dir / "figure_1b_1c_wt_cytokine_contrasts.unfiltered.tsv.gz",
    )
    write_dataframe_gzip(
        wt_combined.loc[
            wt_combined["adjusted_p_value"].lt(0.05)
        ].reset_index(drop=True),
        figure_input_dir / "figure_1b_1c_wt_cytokine_contrasts.fdr05.tsv.gz",
    )
    ko_combined = pd.concat(ko_adapters, ignore_index=True)
    write_dataframe_gzip(
        ko_combined,
        figure_input_dir
        / "figure_2b_s2d_tnfr1_ko1_vs_wt_matched_treatments.unfiltered.tsv.gz",
    )

    manifest_frame = pd.DataFrame(manifest_rows)
    manifest_frame["contrast"] = manifest_frame["contrast"].map(
        lambda value: json.dumps(value, separators=(",", ":"))
    )
    manifest_frame["analysis_status_counts"] = manifest_frame[
        "analysis_status_counts"
    ].map(lambda value: json.dumps(value, sort_keys=True, separators=(",", ":")))
    manifest_frame.to_csv(
        args.output_dir / "analysis_manifest.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )

    runtime = {
        "python": platform.python_version(),
        "pydeseq2": observed_version,
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "n_cpus": args.n_cpus,
        "include_interaction": args.include_interaction,
        "feature_level": "Gene_Symbol after deterministic sum over source Gene_ID",
        "prefilter": "all-zero features within each fitted sample subset only",
        "export_filter": "none; every output contains the full 46,425-symbol universe",
        "output_commit_mode": "staged atomic directory rename",
        "input_sha256": {
            path.name: sha256_file(path)
            for path in (
                args.input_dir / "gene_counts.tsv.gz",
                args.input_dir / "gene_annotations.tsv.gz",
                args.input_dir / "sample_metadata.tsv",
            )
        },
    }
    (args.output_dir / "run_metadata.json").write_text(
        json.dumps(runtime, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    environment = sorted(
        {
            f"{distribution.metadata['Name']}=={distribution.version}"
            for distribution in distributions()
            if distribution.metadata.get("Name")
        },
        key=str.casefold,
    )
    (args.output_dir / "environment.freeze.txt").write_text(
        "\n".join(environment) + "\n", encoding="utf-8"
    )
    staging_output_dir.replace(requested_output_dir)
    atexit.unregister(cleanup_staging)


if __name__ == "__main__":
    main()
