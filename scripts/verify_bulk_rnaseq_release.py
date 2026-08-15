#!/usr/bin/env python3
"""Verify a regenerated bulk RNA-seq run against the validated release.

Floating-point model outputs are compared semantically instead of as compressed
bytes.  Ordinary numeric fields retain at least six significant decimal digits;
p-values are compared in -log10 space so the check remains meaningful across
their full dynamic range.  These bounds allow small BLAS/libm variation across
platforms without accepting scientifically material model drift.
Identifiers, annotations, integer counts, missingness, column order, and all
figure/inference threshold decisions remain exact.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


RELATIVE_TOLERANCE = 1e-6
ABSOLUTE_TOLERANCE = 1e-8
NEG_LOG10_P_TOLERANCE = 1e-4
EXPECTED_RESULT_ROWS = 46_425
EXPECTED_ANALYSES = 10
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
P_VALUE_COLUMNS = frozenset(("p_value", "adjusted_p_value", "pvalue", "padj"))
P_VALUE_FLOOR = 1e-300


class ReleaseVerificationError(ValueError):
    """Raised when a regenerated release violates the validation contract."""


@dataclass(frozen=True)
class DecisionRule:
    """A scientifically meaningful threshold whose membership must be exact."""

    column: str
    operator: str
    threshold: float

    def evaluate(self, values: np.ndarray) -> np.ndarray:
        if self.operator == "lt":
            return values < self.threshold
        if self.operator == "gt":
            return values > self.threshold
        if self.operator == "ge":
            return values >= self.threshold
        raise AssertionError(f"Unsupported decision-rule operator: {self.operator}")

    @property
    def label(self) -> str:
        symbol = {"lt": "<", "gt": ">", "ge": ">="}[self.operator]
        return f"{self.column} {symbol} {self.threshold:g}"


@dataclass(frozen=True)
class TableContract:
    """Exact and floating-point columns for one versioned release table."""

    generated_relative_path: str
    release_filename: str
    key_columns: tuple[str, ...]
    exact_columns: tuple[str, ...]
    numeric_columns: tuple[str, ...]
    decision_rules: tuple[DecisionRule, ...]

    @property
    def columns(self) -> tuple[str, ...]:
        return self.key_columns + self.exact_columns + self.numeric_columns


FIGURE_1_EFFECT = "log2_fold_change_treatment_vs_untreated"
FIGURE_2_EFFECT = "log2_fold_change_ko1_vs_wt"
FIGURE_NUMERIC_PREFIX = ("base_mean",)
FIGURE_NUMERIC_SUFFIX = ("lfc_se", "p_value", "adjusted_p_value")


def figure_decision_rules(effect_column: str) -> tuple[DecisionRule, ...]:
    return (
        DecisionRule("base_mean", "ge", 30.0),
        DecisionRule(effect_column, "gt", 1.0),
        DecisionRule(effect_column, "lt", -1.0),
        DecisionRule("p_value", "lt", 0.05),
        DecisionRule("adjusted_p_value", "lt", 0.05),
    )


TABLE_CONTRACTS = (
    TableContract(
        generated_relative_path=(
            "figure_inputs/"
            "figure_1b_1c_wt_cytokine_contrasts.unfiltered.tsv.gz"
        ),
        release_filename="figure_1b_1c_wt_cytokine_contrasts.unfiltered.tsv.gz",
        key_columns=("condition", "gene_symbol"),
        exact_columns=(),
        numeric_columns=(
            *FIGURE_NUMERIC_PREFIX,
            FIGURE_1_EFFECT,
            "lfc_se",
            "wald_statistic_treatment_vs_untreated",
            *FIGURE_NUMERIC_SUFFIX[1:],
        ),
        decision_rules=figure_decision_rules(FIGURE_1_EFFECT),
    ),
    TableContract(
        generated_relative_path=(
            "figure_inputs/"
            "figure_2b_s2d_tnfr1_ko1_vs_wt_matched_treatments.unfiltered.tsv.gz"
        ),
        release_filename=(
            "figure_2b_s2d_tnfr1_ko1_vs_wt_matched_treatments.unfiltered.tsv.gz"
        ),
        key_columns=("condition", "gene_symbol"),
        exact_columns=(),
        numeric_columns=(
            *FIGURE_NUMERIC_PREFIX,
            FIGURE_2_EFFECT,
            "lfc_se",
            "wald_statistic_ko1_vs_wt",
            *FIGURE_NUMERIC_SUFFIX[1:],
        ),
        decision_rules=figure_decision_rules(FIGURE_2_EFFECT),
    ),
    *(
        TableContract(
            generated_relative_path=filename,
            release_filename=filename,
            key_columns=("Gene_Symbol",),
            exact_columns=(
                "source_gene_ids",
                "n_source_gene_ids",
                "total_count_in_model_subset",
                "analysis_status",
            ),
            numeric_columns=(
                "baseMean",
                "log2FoldChange",
                "lfcSE",
                "stat",
                "pvalue",
                "padj",
            ),
            decision_rules=(
                DecisionRule("baseMean", "ge", 30.0),
                DecisionRule("log2FoldChange", "gt", 1.0),
                DecisionRule("log2FoldChange", "lt", -1.0),
                DecisionRule("pvalue", "lt", 0.05),
                DecisionRule("padj", "lt", 0.05),
            ),
        )
        for filename in (
            "interaction_TNFR1_KO1_vs_WT_IFNG_vs_control.unfiltered.tsv.gz",
            "interaction_TNFR1_KO1_vs_WT_TNF_vs_control.unfiltered.tsv.gz",
            "interaction_TNFR1_KO1_vs_WT_TNF_IFNG_vs_control.unfiltered.tsv.gz",
        )
    ),
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_gzip_payload(path: Path) -> str:
    digest = hashlib.sha256()
    with gzip.open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_text_table(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise ReleaseVerificationError(f"Required release table is missing: {path}")
    return pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        na_filter=False,
    )


def require_columns(frame: pd.DataFrame, expected: Iterable[str], path: Path) -> None:
    expected_columns = list(expected)
    if frame.columns.tolist() != expected_columns:
        raise ReleaseVerificationError(
            f"Column contract differs for {path}: expected {expected_columns!r}, "
            f"observed {frame.columns.tolist()!r}"
        )


def first_exact_difference(
    generated: pd.DataFrame,
    release: pd.DataFrame,
    columns: tuple[str, ...],
) -> tuple[int, str] | None:
    for column in columns:
        different = generated[column].ne(release[column]).to_numpy()
        if different.any():
            return int(np.flatnonzero(different)[0]), column
    return None


def parse_numeric_column(frame: pd.DataFrame, column: str, path: Path) -> np.ndarray:
    tokens = frame[column]
    missing = tokens.eq("NA")
    try:
        parsed = pd.to_numeric(tokens.mask(missing), errors="raise").to_numpy(float)
    except (TypeError, ValueError) as error:
        raise ReleaseVerificationError(
            f"{path}: {column} contains a non-numeric value other than the NA sentinel"
        ) from error
    unexpected_missing = np.isnan(parsed) & ~missing.to_numpy()
    if unexpected_missing.any():
        row = int(np.flatnonzero(unexpected_missing)[0])
        raise ReleaseVerificationError(
            f"{path}: {column} row {row} uses an empty/non-canonical missing value"
        )
    nonmissing = ~missing.to_numpy()
    if not np.isfinite(parsed[nonmissing]).all():
        row = int(np.flatnonzero(nonmissing & ~np.isfinite(parsed))[0])
        raise ReleaseVerificationError(
            f"{path}: {column} row {row} contains a non-finite numeric value"
        )
    return parsed


def compare_numeric_column(
    generated: pd.DataFrame,
    release: pd.DataFrame,
    column: str,
    generated_path: Path,
    release_path: Path,
) -> tuple[np.ndarray, np.ndarray, float, float, float]:
    observed = parse_numeric_column(generated, column, generated_path)
    expected = parse_numeric_column(release, column, release_path)
    observed_missing = np.isnan(observed)
    expected_missing = np.isnan(expected)
    if not np.array_equal(observed_missing, expected_missing):
        row = int(np.flatnonzero(observed_missing != expected_missing)[0])
        raise ReleaseVerificationError(
            f"Missingness differs for {generated_path}, column {column}, row {row}"
        )

    present = ~expected_missing
    observed_present = observed[present]
    expected_present = expected[present]
    absolute_error = np.abs(observed_present - expected_present)
    if column in P_VALUE_COLUMNS:
        if (
            (observed_present < 0).any()
            or (observed_present > 1).any()
            or (expected_present < 0).any()
            or (expected_present > 1).any()
        ):
            raise ReleaseVerificationError(
                f"{generated_path}: {column} contains a value outside [0, 1]"
            )
        observed_scale = -np.log10(
            np.clip(observed_present, P_VALUE_FLOOR, 1.0)
        )
        expected_scale = -np.log10(
            np.clip(expected_present, P_VALUE_FLOOR, 1.0)
        )
        comparison_error = np.abs(observed_scale - expected_scale)
        allowance = np.full_like(comparison_error, NEG_LOG10_P_TOLERANCE)
        tolerance_description = (
            f"max |delta -log10(p)|={NEG_LOG10_P_TOLERANCE:g}"
        )
    else:
        comparison_error = absolute_error
        allowance = (
            ABSOLUTE_TOLERANCE
            + RELATIVE_TOLERANCE * np.abs(expected_present)
        )
        tolerance_description = (
            f"rtol={RELATIVE_TOLERANCE:g}, atol={ABSOLUTE_TOLERANCE:g}"
        )
    outside = comparison_error > allowance
    if outside.any():
        present_rows = np.flatnonzero(present)
        local_row = int(np.flatnonzero(outside)[0])
        row = int(present_rows[local_row])
        raise ReleaseVerificationError(
            f"Numeric drift exceeds {tolerance_description} for {generated_path}, "
            f"column {column}, "
            f"row {row}: regenerated={observed[row]:.17g}, "
            f"release={expected[row]:.17g}, "
            f"comparison_error={comparison_error[local_row]:.3g}, "
            f"allowance={allowance[local_row]:.3g}"
        )

    max_absolute_error = float(absolute_error.max(initial=0.0))
    normalized_error = np.divide(
        comparison_error,
        allowance,
        out=np.zeros_like(comparison_error),
        where=allowance > 0,
    )
    max_normalized_error = float(normalized_error.max(initial=0.0))
    max_neg_log10_p_error = (
        float(comparison_error.max(initial=0.0))
        if column in P_VALUE_COLUMNS
        else 0.0
    )
    return (
        observed,
        expected,
        max_absolute_error,
        max_normalized_error,
        max_neg_log10_p_error,
    )


def compare_table(
    generated_path: Path,
    release_path: Path,
    contract: TableContract,
) -> dict[str, object]:
    generated = load_text_table(generated_path)
    release = load_text_table(release_path)
    generated_compressed_sha256 = sha256_file(generated_path)
    release_compressed_sha256 = sha256_file(release_path)
    generated_payload_sha256 = sha256_gzip_payload(generated_path)
    release_payload_sha256 = sha256_gzip_payload(release_path)
    require_columns(generated, contract.columns, generated_path)
    require_columns(release, contract.columns, release_path)
    if len(generated) != len(release):
        raise ReleaseVerificationError(
            f"Row count differs for {generated_path}: regenerated={len(generated)}, "
            f"release={len(release)}"
        )

    keys = list(contract.key_columns)
    for frame, path in ((generated, generated_path), (release, release_path)):
        if frame[keys].eq("").any(axis=None):
            raise ReleaseVerificationError(f"Semantic key is empty in {path}")
        duplicate = frame.duplicated(keys, keep=False)
        if duplicate.any():
            row = int(np.flatnonzero(duplicate.to_numpy())[0])
            raise ReleaseVerificationError(
                f"Semantic key is duplicated in {path} at row {row}"
            )

    exact_columns = contract.key_columns + contract.exact_columns
    difference = first_exact_difference(generated, release, exact_columns)
    if difference is not None:
        row, column = difference
        key = generated.loc[row, keys].to_dict()
        raise ReleaseVerificationError(
            f"Exact semantic field differs for {generated_path}, key={key!r}, "
            f"column={column}: regenerated={generated.at[row, column]!r}, "
            f"release={release.at[row, column]!r}"
        )

    parsed: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    max_absolute_error = 0.0
    max_normalized_error = 0.0
    max_neg_log10_p_error = 0.0
    for column in contract.numeric_columns:
        (
            observed,
            expected,
            column_absolute,
            column_normalized,
            column_neg_log10_p,
        ) = compare_numeric_column(
            generated, release, column, generated_path, release_path
        )
        parsed[column] = (observed, expected)
        max_absolute_error = max(max_absolute_error, column_absolute)
        max_normalized_error = max(max_normalized_error, column_normalized)
        max_neg_log10_p_error = max(
            max_neg_log10_p_error, column_neg_log10_p
        )

    for rule in contract.decision_rules:
        observed, expected = parsed[rule.column]
        observed_decision = rule.evaluate(observed)
        expected_decision = rule.evaluate(expected)
        if not np.array_equal(observed_decision, expected_decision):
            row = int(np.flatnonzero(observed_decision != expected_decision)[0])
            key = generated.loc[row, keys].to_dict()
            raise ReleaseVerificationError(
                f"Scientific decision differs for {generated_path}, key={key!r}, "
                f"rule={rule.label}: regenerated={observed[row]:.17g}, "
                f"release={expected[row]:.17g}"
            )

    return {
        "file": contract.release_filename,
        "rows": len(generated),
        "generated_bytes": generated_path.stat().st_size,
        "release_bytes": release_path.stat().st_size,
        "generated_compressed_sha256": generated_compressed_sha256,
        "release_compressed_sha256": release_compressed_sha256,
        "generated_decompressed_sha256": generated_payload_sha256,
        "release_decompressed_sha256": release_payload_sha256,
        "compression_only_difference": (
            generated_compressed_sha256 != release_compressed_sha256
            and generated_payload_sha256 == release_payload_sha256
        ),
        "exact_schema_keys_order_text_and_na_masks": True,
        "decision_invariants": [rule.label for rule in contract.decision_rules],
        "max_absolute_error": max_absolute_error,
        "max_neg_log10_p_error": max_neg_log10_p_error,
        "max_tolerance_fraction": max_normalized_error,
    }


def compare_manifest(generated_path: Path, release_path: Path) -> pd.DataFrame:
    generated = load_text_table(generated_path)
    release = load_text_table(release_path)
    require_columns(generated, release.columns, generated_path)
    if "sha256" not in generated.columns:
        raise ReleaseVerificationError("analysis_manifest.tsv lacks the sha256 column")
    for frame, path in ((generated, generated_path), (release, release_path)):
        invalid_hash = ~frame["sha256"].map(
            lambda value: SHA256_PATTERN.fullmatch(value) is not None
        )
        if invalid_hash.any():
            row = int(np.flatnonzero(invalid_hash.to_numpy())[0])
            raise ReleaseVerificationError(
                f"Invalid output SHA-256 in {path} at row {row}"
            )

    comparison_columns = tuple(
        column for column in release.columns if column != "sha256"
    )
    if len(generated) != len(release):
        raise ReleaseVerificationError(
            "analysis_manifest.tsv row count differs from the validated release"
        )
    difference = first_exact_difference(generated, release, comparison_columns)
    if difference is not None:
        row, column = difference
        raise ReleaseVerificationError(
            f"analysis_manifest.tsv differs outside output SHA hashes at row {row}, "
            f"column {column}: regenerated={generated.at[row, column]!r}, "
            f"release={release.at[row, column]!r}"
        )
    return generated


def manifest_hash_provenance(
    generated_path: Path, release_path: Path
) -> list[dict[str, object]]:
    generated = load_text_table(generated_path)
    release = load_text_table(release_path)
    return [
        {
            "analysis_id": generated.at[row, "analysis_id"],
            "output_file": generated.at[row, "output_file"],
            "generated_sha256": generated.at[row, "sha256"],
            "release_sha256": release.at[row, "sha256"],
            "hashes_equal": (
                generated.at[row, "sha256"] == release.at[row, "sha256"]
            ),
        }
        for row in range(len(generated))
    ]


def verify_complete_exports(generated_dir: Path, manifest: pd.DataFrame) -> None:
    if len(manifest) != EXPECTED_ANALYSES or manifest["analysis_id"].duplicated().any():
        raise ReleaseVerificationError(
            f"Expected {EXPECTED_ANALYSES} unique primary/interaction analyses"
        )
    for row in manifest.itertuples(index=False):
        relative_path = Path(row.output_file)
        if relative_path.is_absolute() or ".." in relative_path.parts:
            raise ReleaseVerificationError(
                f"Unsafe output path in analysis_manifest.tsv: {row.output_file!r}"
            )
        path = generated_dir / relative_path
        if not path.is_file():
            raise ReleaseVerificationError(f"Manifest output is missing: {path}")
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            header = handle.readline().rstrip("\n").split("\t")
            observed_rows = sum(1 for _ in handle)
        if observed_rows != EXPECTED_RESULT_ROWS or "Gene_Symbol" not in header:
            raise ReleaseVerificationError(
                f"Unexpected unfiltered result contract: {path}"
            )
        if int(row.rows) != observed_rows:
            raise ReleaseVerificationError(
                f"Manifest row count does not match {path}: "
                f"manifest={row.rows}, observed={observed_rows}"
            )
        observed_sha256 = sha256_file(path)
        if observed_sha256 != row.sha256:
            raise ReleaseVerificationError(
                f"Manifest output SHA-256 does not match {path}"
            )


def compare_metadata(generated_path: Path, release_path: Path) -> None:
    try:
        generated = json.loads(generated_path.read_text(encoding="utf-8"))
        release = json.loads(release_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ReleaseVerificationError(f"Cannot read run metadata: {error}") from error
    if generated != release:
        differing_keys = sorted(
            key
            for key in set(generated) | set(release)
            if generated.get(key) != release.get(key)
        )
        raise ReleaseVerificationError(
            f"run_metadata.json differs from the validated release: {differing_keys!r}"
        )


def compare_environment(generated_path: Path, release_path: Path) -> None:
    try:
        generated = generated_path.read_text(encoding="utf-8")
        release = release_path.read_text(encoding="utf-8")
    except OSError as error:
        raise ReleaseVerificationError(f"Cannot read environment freeze: {error}") from error
    if generated != release:
        raise ReleaseVerificationError(
            "environment.freeze.txt differs from the validated release"
        )


def verify_release(generated_dir: Path, release_dir: Path) -> list[dict[str, object]]:
    manifest = compare_manifest(
        generated_dir / "analysis_manifest.tsv",
        release_dir / "analysis_manifest.tsv",
    )
    verify_complete_exports(generated_dir, manifest)
    compare_metadata(
        generated_dir / "run_metadata.json",
        release_dir / "run_metadata.json",
    )
    compare_environment(
        generated_dir / "environment.freeze.txt",
        release_dir / "environment.freeze.txt",
    )
    summaries = []
    for contract in TABLE_CONTRACTS:
        summaries.append(
            compare_table(
                generated_dir / contract.generated_relative_path,
                release_dir / contract.release_filename,
                contract,
            )
        )
    return summaries


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generated-dir", type=Path, required=True)
    parser.add_argument("--release-dir", type=Path, required=True)
    parser.add_argument(
        "--report",
        type=Path,
        help="Write a JSON comparison report, including failures, for CI artifacts.",
    )
    args = parser.parse_args()
    report: dict[str, object] = {
        "schema_version": 1,
        "generated_dir": str(args.generated_dir),
        "release_dir": str(args.release_dir),
        "tolerances": {
            "ordinary_numeric_rtol": RELATIVE_TOLERANCE,
            "ordinary_numeric_atol": ABSOLUTE_TOLERANCE,
            "max_absolute_delta_neg_log10_p": NEG_LOG10_P_TOLERANCE,
        },
        "status": "failed",
        "exact_checks": {
            "manifest_excluding_output_sha256": "not_completed",
            "generated_output_sha256_self_consistency": "not_completed",
            "run_metadata": "not_completed",
            "environment_freeze": "not_completed",
        },
        "tables": [],
    }
    error: Exception | None = None
    try:
        summaries = verify_release(args.generated_dir, args.release_dir)
        report["status"] = "passed"
        report["exact_checks"] = {
            "manifest_excluding_output_sha256": "passed",
            "generated_output_sha256_self_consistency": "passed",
            "run_metadata": "passed",
            "environment_freeze": "passed",
        }
        report["tables"] = summaries
        report["manifest_output_hash_provenance"] = manifest_hash_provenance(
            args.generated_dir / "analysis_manifest.tsv",
            args.release_dir / "analysis_manifest.tsv",
        )
    except Exception as caught:
        error = caught
        report["error_type"] = type(caught).__name__
        report["error"] = str(caught)
    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    if error is not None:
        raise SystemExit(str(error)) from error
    for summary in summaries:
        print(
            "SEMANTIC_MATCH "
            f"{summary['file']}: rows={summary['rows']}, exact_contract=ok, "
            f"max_abs_error={summary['max_absolute_error']:.3g}, "
            f"max_delta_neg_log10_p={summary['max_neg_log10_p_error']:.3g}, "
            f"max_tolerance_fraction={summary['max_tolerance_fraction']:.3g}"
        )
    print(
        "BULK_RNASEQ_RESULTS_OK "
        f"(rtol={RELATIVE_TOLERANCE:g}, atol={ABSOLUTE_TOLERANCE:g}, "
        f"max_delta_neg_log10_p={NEG_LOG10_P_TOLERANCE:g})"
    )


if __name__ == "__main__":
    main()
