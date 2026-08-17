#!/usr/bin/env python3
"""Validate aggregate TCR projection outputs without exposing cell identifiers."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


EXPECTED_CLUSTERS = ("C0", "C1", "C2", "C3", "C4", "C5")
EXPECTED_LABELS = {
    "C0": "Mixed CD4/KLRB1-associated activated state",
    "C1": "Cycling T-cell state I",
    "C2": "Cytokine-expressing effector state",
    "C3": "CD8/TRDC-associated cytotoxic state",
    "C4": "Cycling T-cell state II",
    "C5": "CCR7/IL7R/HLA-II-associated state",
}
EXPECTED_CYCLE_GENES = frozenset(
    {"TK1", "MKI67", "AURKB", "TOP2A", "UBE2C", "HMGB2", "TYMS", "HMMR", "PTTG2"}
)
BY_CLUSTER_COLUMNS = (
    "cluster_short",
    "cluster_annotation",
    "n_cells",
    "c6_full_rank_score_mean",
    "c6_full_rank_score_median",
    "c6_full_rank_score_q25",
    "c6_full_rank_score_q75",
    "c6_cycle_rank_score_mean",
    "c6_cycle_rank_score_median",
    "c6_cycle_rank_score_q25",
    "c6_cycle_rank_score_q75",
    "c6_noncycle_rank_score_mean",
    "c6_noncycle_rank_score_median",
    "c6_noncycle_rank_score_q25",
    "c6_noncycle_rank_score_q75",
    "cxcl13_detected_cells",
    "cxcl13_detection_fraction",
    "cxcl13_mean_log_normalized_expression",
)
GENE_COVERAGE_COLUMNS = (
    "gene",
    "component",
    "frozen_avg_log2FC",
    "in_targeted_panel",
    "used_in_score",
    "exclusion_reason",
)
CLUSTER_COUNT_COLUMNS = ("cluster_short", "n_cells", "fraction")


class ContractError(ValueError):
    """Raised when an output violates the release contract."""


def read_delimited(path: Path, delimiter: str) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise ContractError(f"Missing or empty file: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ContractError(f"Missing header: {path}")
        rows = list(reader)
    return reader.fieldnames, rows


def require_columns(path: Path, observed: list[str], expected: tuple[str, ...]) -> None:
    if tuple(observed) != expected:
        raise ContractError(
            f"Unexpected columns in {path}: {observed}; expected {list(expected)}"
        )


def parse_int(value: str, *, field: str, row_name: str) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError) as error:
        raise ContractError(f"{row_name}: {field} must be an integer") from error
    if parsed < 0:
        raise ContractError(f"{row_name}: {field} must be non-negative")
    return parsed


def parse_float(value: str, *, field: str, row_name: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as error:
        raise ContractError(f"{row_name}: {field} must be numeric") from error
    if not math.isfinite(parsed):
        raise ContractError(f"{row_name}: {field} must be finite")
    return parsed


def parse_bool(value: str, *, field: str, row_name: str) -> bool:
    normalized = value.strip().lower()
    if normalized == "true":
        return True
    if normalized == "false":
        return False
    raise ContractError(f"{row_name}: {field} must be TRUE or FALSE")


def read_c6_signature(path: Path) -> list[tuple[str, float]]:
    header, rows = read_delimited(path, ",")
    required = {"cluster", "gene", "avg_log2FC"}
    if not required.issubset(header):
        raise ContractError(f"Signature resource lacks {sorted(required)}: {path}")
    genes = [
        (
            row["gene"].strip(),
            parse_float(
                row["avg_log2FC"],
                field="avg_log2FC",
                row_name=row["gene"].strip(),
            ),
        )
        for row in rows
        if row["cluster"].strip() == "6"
    ]
    if len(genes) != 20 or len({gene for gene, _ in genes}) != 20:
        raise ContractError("Frozen tumor-co-culture C6 signature must contain 20 genes")
    return genes


def read_cluster_counts(path: Path) -> dict[str, int]:
    header, rows = read_delimited(path, ",")
    require_columns(path, header, CLUSTER_COUNT_COLUMNS)
    observed: dict[str, int] = {}
    for row in rows:
        cluster = row["cluster_short"].strip()
        if cluster in observed:
            raise ContractError(f"Duplicate cluster in {path}: {cluster}")
        observed[cluster] = parse_int(
            row["n_cells"], field="n_cells", row_name=cluster
        )
        fraction = parse_float(row["fraction"], field="fraction", row_name=cluster)
        if not 0.0 <= fraction <= 1.0:
            raise ContractError(f"{cluster}: fraction must be in [0, 1]")
    if tuple(sorted(observed)) != EXPECTED_CLUSTERS:
        raise ContractError(f"Cluster-count table must contain C0-C5: {path}")
    total = sum(observed.values())
    if total <= 0:
        raise ContractError("Cluster-count table has no retained cells")
    for row in rows:
        cluster = row["cluster_short"].strip()
        expected_fraction = observed[cluster] / total
        observed_fraction = float(row["fraction"])
        if not math.isclose(observed_fraction, expected_fraction, rel_tol=1e-10, abs_tol=1e-12):
            raise ContractError(f"{cluster}: fraction does not match n_cells / total")
    return observed


def verify_by_cluster(path: Path, expected_counts: dict[str, int]) -> None:
    header, rows = read_delimited(path, "\t")
    require_columns(path, header, BY_CLUSTER_COLUMNS)
    if len(rows) != len(EXPECTED_CLUSTERS):
        raise ContractError(f"Projection summary must contain six cluster rows: {path}")
    seen: set[str] = set()
    for row in rows:
        cluster = row["cluster_short"].strip()
        if cluster in seen:
            raise ContractError(f"Duplicate projection cluster: {cluster}")
        seen.add(cluster)
        if cluster not in EXPECTED_LABELS:
            raise ContractError(f"Unexpected projection cluster: {cluster}")
        if row["cluster_annotation"].strip() != EXPECTED_LABELS[cluster]:
            raise ContractError(f"{cluster}: annotation differs from the frozen TCR label")
        n_cells = parse_int(row["n_cells"], field="n_cells", row_name=cluster)
        if n_cells != expected_counts[cluster]:
            raise ContractError(f"{cluster}: projection n_cells differs from count table")

        for score_prefix in (
            "c6_full_rank_score",
            "c6_cycle_rank_score",
            "c6_noncycle_rank_score",
        ):
            q25 = parse_float(
                row[f"{score_prefix}_q25"],
                field=f"{score_prefix}_q25",
                row_name=cluster,
            )
            median = parse_float(
                row[f"{score_prefix}_median"],
                field=f"{score_prefix}_median",
                row_name=cluster,
            )
            q75 = parse_float(
                row[f"{score_prefix}_q75"],
                field=f"{score_prefix}_q75",
                row_name=cluster,
            )
            mean = parse_float(
                row[f"{score_prefix}_mean"],
                field=f"{score_prefix}_mean",
                row_name=cluster,
            )
            if not q25 <= median <= q75:
                raise ContractError(
                    f"{cluster}: {score_prefix} quantiles are not ordered"
                )
            if not all(0.0 <= value <= 1.0 for value in (mean, q25, median, q75)):
                raise ContractError(
                    f"{cluster}: {score_prefix} summaries must be rank scores in [0, 1]"
                )

        detected = parse_int(
            row["cxcl13_detected_cells"],
            field="cxcl13_detected_cells",
            row_name=cluster,
        )
        if detected > n_cells:
            raise ContractError(f"{cluster}: CXCL13-detected cells exceed n_cells")
        detection_fraction = parse_float(
            row["cxcl13_detection_fraction"],
            field="cxcl13_detection_fraction",
            row_name=cluster,
        )
        if not 0.0 <= detection_fraction <= 1.0:
            raise ContractError(f"{cluster}: CXCL13 detection fraction must be in [0, 1]")
        if not math.isclose(
            detection_fraction,
            detected / n_cells,
            rel_tol=1e-10,
            abs_tol=1e-12,
        ):
            raise ContractError(
                f"{cluster}: CXCL13 detection fraction differs from detected/n_cells"
            )
        mean_expression = parse_float(
            row["cxcl13_mean_log_normalized_expression"],
            field="cxcl13_mean_log_normalized_expression",
            row_name=cluster,
        )
        if mean_expression < 0:
            raise ContractError(
                f"{cluster}: mean log-normalized CXCL13 expression must be non-negative"
            )
    if tuple(sorted(seen)) != EXPECTED_CLUSTERS:
        raise ContractError("Projection summary must contain exactly C0-C5")


def verify_gene_coverage(
    path: Path, expected_signature: list[tuple[str, float]]
) -> None:
    header, rows = read_delimited(path, "\t")
    require_columns(path, header, GENE_COVERAGE_COLUMNS)
    if len(rows) != len(expected_signature):
        raise ContractError("Gene-coverage table must have one row per frozen C6 gene")
    expected_genes = [gene for gene, _ in expected_signature]
    observed_genes = [row["gene"].strip() for row in rows]
    if observed_genes != expected_genes:
        raise ContractError(
            "Gene-coverage rows must preserve the frozen C6 signature order"
        )
    if len(set(observed_genes)) != len(observed_genes):
        raise ContractError("Gene-coverage table contains duplicate genes")
    expected_fc = dict(expected_signature)
    for row in rows:
        gene = row["gene"].strip()
        component = row["component"].strip()
        expected_component = "cycle" if gene in EXPECTED_CYCLE_GENES else "noncycle"
        if component != expected_component:
            raise ContractError(
                f"{gene}: component must be {expected_component}, observed {component}"
            )
        frozen_fc = parse_float(
            row["frozen_avg_log2FC"], field="frozen_avg_log2FC", row_name=gene
        )
        if not math.isclose(
            frozen_fc, expected_fc[gene], rel_tol=1e-12, abs_tol=1e-12
        ):
            raise ContractError(
                f"{gene}: frozen_avg_log2FC differs from the frozen signature"
            )
        in_panel = parse_bool(
            row["in_targeted_panel"], field="in_targeted_panel", row_name=gene
        )
        used = parse_bool(row["used_in_score"], field="used_in_score", row_name=gene)
        reason = row["exclusion_reason"].strip()
        if not in_panel or not used:
            raise ContractError(
                f"{gene}: every frozen C6 gene is present and must be used in this panel"
            )
        if reason:
            raise ContractError(f"{gene}: exclusion_reason must be blank when used")
    cxcl13 = rows[expected_genes.index("CXCL13")]
    if cxcl13["component"].strip() != "noncycle":
        raise ContractError("CXCL13 must remain in the noncycle component")
    if not parse_bool(
        cxcl13["used_in_score"], field="used_in_score", row_name="CXCL13"
    ):
        raise ContractError("CXCL13 must be included in the C6 projection score")


def verify_release(
    by_cluster: Path,
    gene_coverage: Path,
    cluster_counts: Path,
    signature_resource: Path,
) -> None:
    expected_counts = read_cluster_counts(cluster_counts)
    expected_signature = read_c6_signature(signature_resource)
    verify_by_cluster(by_cluster, expected_counts)
    verify_gene_coverage(gene_coverage, expected_signature)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--by-cluster", type=Path, required=True)
    parser.add_argument("--gene-coverage", type=Path, required=True)
    parser.add_argument("--cluster-counts", type=Path, required=True)
    parser.add_argument("--signature-resource", type=Path, required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        verify_release(
            args.by_cluster,
            args.gene_coverage,
            args.cluster_counts,
            args.signature_resource,
        )
    except ContractError as error:
        raise SystemExit(f"TCR C6 projection contract failed: {error}") from error
    print("TCR C6 projection contract passed (aggregate outputs only).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
