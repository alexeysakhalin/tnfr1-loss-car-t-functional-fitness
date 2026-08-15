#!/usr/bin/env python3
"""Verify regenerated bulk RNA-seq results against the validated release.

The fitted values produced by DESeq2-family software can vary in numerically
unstable tails when the same locked environment is run with a different BLAS or
libm implementation. This verifier treats raw numeric differences as
diagnostics, not as a global pass/fail tolerance. It keeps the parts that define
the analysis and manuscript conclusions exact: schema, row keys and order,
annotations, integer fields, missingness, analysis metadata, threshold
membership, combined DEG calls, Venn membership, and prespecified interaction
gene conclusions. Fold changes for genes displayed as figure labels have a
narrow manuscript-precision guard. Each table must also satisfy its own
Wald = log2FC / SE identity.
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


EXPECTED_RESULT_ROWS = 46_425
EXPECTED_ANALYSES = 10
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
INTEGER_PATTERN = re.compile(r"0|[1-9][0-9]*")
P_VALUE_COLUMNS = frozenset(("p_value", "adjusted_p_value", "pvalue", "padj"))
P_VALUE_FLOOR = 1e-300
WALD_RTOL = 1e-10
WALD_ATOL = 1e-12
PRESPECIFIED_INTERACTION_GENES = ("ICAM1", "MLKL", "GSDME", "IRF1")
MANUSCRIPT_LABEL_LFC_ATOL = 0.001
FIGURE_1_LABEL_GENES = (
    "ICAM1",
    "IRF1",
    "AIM2",
    "CASP1",
    "CASP4",
    "MLKL",
    "BAK1",
    "CASP7",
    "FAS",
    "CASP8",
)
FIGURE_2_LABEL_GENES = ("ICAM1", "MLKL", "GSDME", "IRF1")
MANUSCRIPT_LABEL_GENES = tuple(
    dict.fromkeys(FIGURE_1_LABEL_GENES + FIGURE_2_LABEL_GENES)
)


class ReleaseVerificationError(ValueError):
    """Raised when a regenerated release violates the validation contract."""


@dataclass(frozen=True)
class TableContract:
    """Column and scientific-outcome contract for one release table."""

    generated_relative_path: str
    release_filename: str
    key_columns: tuple[str, ...]
    exact_columns: tuple[str, ...]
    integer_columns: tuple[str, ...]
    numeric_columns: tuple[str, ...]
    base_mean_column: str
    effect_column: str
    standard_error_column: str
    wald_column: str
    adjusted_p_column: str
    venn_conditions: tuple[str, ...] = ()
    venn_direction: str | None = None
    prespecified_genes: tuple[str, ...] = ()
    guarded_effect_genes: tuple[str, ...] = ()

    @property
    def columns(self) -> tuple[str, ...]:
        return self.key_columns + self.exact_columns + self.numeric_columns


FIGURE_1_EFFECT = "log2_fold_change_treatment_vs_untreated"
FIGURE_1_WALD = "wald_statistic_treatment_vs_untreated"
FIGURE_2_EFFECT = "log2_fold_change_ko1_vs_wt"
FIGURE_2_WALD = "wald_statistic_ko1_vs_wt"
CYTOKINE_VENN_CONDITIONS = ("TNF", "IFNG", "TNF_IFNG")


TABLE_CONTRACTS = (
    TableContract(
        generated_relative_path=(
            "figure_inputs/"
            "figure_1b_1c_wt_cytokine_contrasts.unfiltered.tsv.gz"
        ),
        release_filename="figure_1b_1c_wt_cytokine_contrasts.unfiltered.tsv.gz",
        key_columns=("condition", "gene_symbol"),
        exact_columns=(),
        integer_columns=(),
        numeric_columns=(
            "base_mean",
            FIGURE_1_EFFECT,
            "lfc_se",
            FIGURE_1_WALD,
            "p_value",
            "adjusted_p_value",
        ),
        base_mean_column="base_mean",
        effect_column=FIGURE_1_EFFECT,
        standard_error_column="lfc_se",
        wald_column=FIGURE_1_WALD,
        adjusted_p_column="adjusted_p_value",
        venn_conditions=CYTOKINE_VENN_CONDITIONS,
        venn_direction="up",
        guarded_effect_genes=FIGURE_1_LABEL_GENES,
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
        integer_columns=(),
        numeric_columns=(
            "base_mean",
            FIGURE_2_EFFECT,
            "lfc_se",
            FIGURE_2_WALD,
            "p_value",
            "adjusted_p_value",
        ),
        base_mean_column="base_mean",
        effect_column=FIGURE_2_EFFECT,
        standard_error_column="lfc_se",
        wald_column=FIGURE_2_WALD,
        adjusted_p_column="adjusted_p_value",
        venn_conditions=CYTOKINE_VENN_CONDITIONS,
        venn_direction="down",
        guarded_effect_genes=FIGURE_2_LABEL_GENES,
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
            integer_columns=(
                "n_source_gene_ids",
                "total_count_in_model_subset",
            ),
            numeric_columns=(
                "baseMean",
                "log2FoldChange",
                "lfcSE",
                "stat",
                "pvalue",
                "padj",
            ),
            base_mean_column="baseMean",
            effect_column="log2FoldChange",
            standard_error_column="lfcSE",
            wald_column="stat",
            adjusted_p_column="padj",
            prespecified_genes=PRESPECIFIED_INTERACTION_GENES,
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


def sha256_strings(values: Iterable[str]) -> str:
    payload = "\n".join(sorted(values)).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


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


def validate_integer_columns(
    frame: pd.DataFrame,
    columns: tuple[str, ...],
    path: Path,
) -> None:
    for column in columns:
        invalid = ~frame[column].map(
            lambda value: INTEGER_PATTERN.fullmatch(value) is not None
        )
        if invalid.any():
            row = int(np.flatnonzero(invalid.to_numpy())[0])
            raise ReleaseVerificationError(
                f"{path}: {column} row {row} is not a canonical non-negative integer"
            )


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


def validate_numeric_domain(column: str, values: np.ndarray, path: Path) -> None:
    present = np.isfinite(values)
    if column in P_VALUE_COLUMNS and (
        (values[present] < 0).any() or (values[present] > 1).any()
    ):
        raise ReleaseVerificationError(f"{path}: {column} contains a value outside [0, 1]")
    if column in ("base_mean", "baseMean") and (values[present] < 0).any():
        raise ReleaseVerificationError(f"{path}: {column} contains a negative value")
    if column in ("lfc_se", "lfcSE") and (values[present] <= 0).any():
        raise ReleaseVerificationError(
            f"{path}: {column} must be positive wherever it is reported"
        )


def numeric_delta_summary(
    observed: np.ndarray,
    expected: np.ndarray,
    column: str,
    generated: pd.DataFrame,
    contract: TableContract,
) -> dict[str, object]:
    present = np.isfinite(expected)
    present_rows = np.flatnonzero(present)
    observed_present = observed[present]
    expected_present = expected[present]
    absolute = np.abs(observed_present - expected_present)
    scale = np.maximum(
        np.maximum(np.abs(observed_present), np.abs(expected_present)),
        np.finfo(float).tiny,
    )
    relative = absolute / scale
    different = observed_present != expected_present
    if absolute.size:
        maximum_index = int(np.argmax(absolute))
        maximum_row = int(present_rows[maximum_index])
        max_values = {
            "row_index": maximum_row,
            "key": key_at(generated, maximum_row, contract),
            "generated": float(observed_present[maximum_index]),
            "release": float(expected_present[maximum_index]),
        }
    else:
        max_values = {"generated": None, "release": None}
    summary: dict[str, object] = {
        "compared_nonmissing_values": int(absolute.size),
        "different_numeric_values": int(different.sum()),
        "max_absolute_delta": float(absolute.max(initial=0.0)),
        "median_absolute_delta": (
            float(np.median(absolute)) if absolute.size else 0.0
        ),
        "p95_absolute_delta": (
            float(np.quantile(absolute, 0.95)) if absolute.size else 0.0
        ),
        "p99_absolute_delta": (
            float(np.quantile(absolute, 0.99)) if absolute.size else 0.0
        ),
        "max_symmetric_relative_delta": float(relative.max(initial=0.0)),
        "values_at_max_absolute_delta": max_values,
    }
    if column in P_VALUE_COLUMNS:
        observed_log = -np.log10(np.clip(observed_present, P_VALUE_FLOOR, 1.0))
        expected_log = -np.log10(np.clip(expected_present, P_VALUE_FLOOR, 1.0))
        summary["max_absolute_delta_neg_log10"] = float(
            np.abs(observed_log - expected_log).max(initial=0.0)
        )
    return summary


def parse_and_compare_numeric_columns(
    generated: pd.DataFrame,
    release: pd.DataFrame,
    contract: TableContract,
    generated_path: Path,
    release_path: Path,
) -> tuple[
    dict[str, tuple[np.ndarray, np.ndarray]],
    dict[str, dict[str, object]],
]:
    parsed: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    diagnostics: dict[str, dict[str, object]] = {}
    for column in contract.numeric_columns:
        observed = parse_numeric_column(generated, column, generated_path)
        expected = parse_numeric_column(release, column, release_path)
        observed_missing = np.isnan(observed)
        expected_missing = np.isnan(expected)
        if not np.array_equal(observed_missing, expected_missing):
            row = int(np.flatnonzero(observed_missing != expected_missing)[0])
            raise ReleaseVerificationError(
                f"Missingness differs for {generated_path}, column {column}, row {row}"
            )
        validate_numeric_domain(column, observed, generated_path)
        validate_numeric_domain(column, expected, release_path)
        parsed[column] = (observed, expected)
        diagnostics[column] = numeric_delta_summary(
            observed, expected, column, generated, contract
        )
    return parsed, diagnostics


def key_at(frame: pd.DataFrame, row: int, contract: TableContract) -> dict[str, str]:
    return frame.loc[row, list(contract.key_columns)].to_dict()


def require_equal_mask(
    observed: np.ndarray,
    expected: np.ndarray,
    label: str,
    generated: pd.DataFrame,
    contract: TableContract,
    generated_values: np.ndarray,
    release_values: np.ndarray,
    generated_path: Path,
) -> None:
    if np.array_equal(observed, expected):
        return
    row = int(np.flatnonzero(observed != expected)[0])
    raise ReleaseVerificationError(
        f"Scientific outcome differs for {generated_path}, "
        f"key={key_at(generated, row, contract)!r}, rule={label}: "
        f"regenerated={generated_values[row]:.17g}, "
        f"release={release_values[row]:.17g}"
    )


def combined_category(
    base_mean: np.ndarray,
    effect: np.ndarray,
    adjusted_p: np.ndarray,
) -> np.ndarray:
    category = np.zeros(len(base_mean), dtype=np.int8)
    eligible = (base_mean >= 30.0) & (adjusted_p < 0.05)
    category[eligible & (effect > 1.0)] = 1
    category[eligible & (effect < -1.0)] = -1
    return category


def standalone_padj_flips(
    generated: pd.DataFrame,
    contract: TableContract,
    observed_adjusted_p: np.ndarray,
    expected_adjusted_p: np.ndarray,
    observed_category: np.ndarray,
    expected_category: np.ndarray,
) -> dict[str, object]:
    observed = observed_adjusted_p < 0.05
    expected = expected_adjusted_p < 0.05
    rows = np.flatnonzero(observed != expected)
    records = []
    for row in rows:
        records.append(
            {
                "key": key_at(generated, int(row), contract),
                "generated_adjusted_p": float(observed_adjusted_p[row]),
                "release_adjusted_p": float(expected_adjusted_p[row]),
                "generated_combined_category": int(observed_category[row]),
                "release_combined_category": int(expected_category[row]),
            }
        )
    return {
        "threshold": 0.05,
        "count": int(len(rows)),
        "outcome_categories_unchanged": bool(
            np.array_equal(observed_category[rows], expected_category[rows])
        ),
        "rows": records,
    }


def validate_wald_identity(
    frame: pd.DataFrame,
    values: dict[str, np.ndarray],
    contract: TableContract,
    path: Path,
) -> dict[str, object]:
    effect = values[contract.effect_column]
    standard_error = values[contract.standard_error_column]
    wald = values[contract.wald_column]
    presence = tuple(np.isfinite(array) for array in (effect, standard_error, wald))
    if not (
        np.array_equal(presence[0], presence[1])
        and np.array_equal(presence[0], presence[2])
    ):
        mismatch = (presence[0] != presence[1]) | (presence[0] != presence[2])
        row = int(np.flatnonzero(mismatch)[0])
        raise ReleaseVerificationError(
            f"Incomplete Wald/log2FC/SE triplet in {path}, "
            f"key={key_at(frame, row, contract)!r}"
        )
    complete = presence[0]
    expected_wald = effect[complete] / standard_error[complete]
    observed_wald = wald[complete]
    close = np.isclose(
        observed_wald,
        expected_wald,
        rtol=WALD_RTOL,
        atol=WALD_ATOL,
    )
    if not close.all():
        complete_rows = np.flatnonzero(complete)
        local_row = int(np.flatnonzero(~close)[0])
        row = int(complete_rows[local_row])
        raise ReleaseVerificationError(
            f"Wald statistic is not log2FC/SE in {path}, "
            f"key={key_at(frame, row, contract)!r}: "
            f"wald={wald[row]:.17g}, log2FC/SE={effect[row] / standard_error[row]:.17g}"
        )
    error = np.abs(observed_wald - expected_wald)
    return {
        "complete_triplets": int(complete.sum()),
        "rtol": WALD_RTOL,
        "atol": WALD_ATOL,
        "max_absolute_identity_error": float(error.max(initial=0.0)),
    }


def venn_membership_summary(
    frame: pd.DataFrame,
    category: np.ndarray,
    contract: TableContract,
) -> tuple[dict[str, set[str]], dict[str, object]]:
    direction_value = 1 if contract.venn_direction == "up" else -1
    sets: dict[str, set[str]] = {}
    for condition in contract.venn_conditions:
        selected = (
            frame["condition"].eq(condition).to_numpy()
            & (category == direction_value)
        )
        sets[condition] = set(frame.loc[selected, "gene_symbol"])

    union = set().union(*sets.values())
    regions: dict[str, dict[str, object]] = {}
    for bits in range(1, 1 << len(contract.venn_conditions)):
        included = tuple(
            condition
            for index, condition in enumerate(contract.venn_conditions)
            if bits & (1 << index)
        )
        excluded = set(contract.venn_conditions).difference(included)
        members = {
            gene
            for gene in union
            if all(gene in sets[condition] for condition in included)
            and all(gene not in sets[condition] for condition in excluded)
        }
        label = " & ".join(included) + " only"
        regions[label] = {
            "count": len(members),
            "membership_sha256": sha256_strings(members),
        }
    summary = {
        "conditions": list(contract.venn_conditions),
        "direction": contract.venn_direction,
        "sets": {
            condition: {
                "count": len(sets[condition]),
                "membership_sha256": sha256_strings(sets[condition]),
            }
            for condition in contract.venn_conditions
        },
        "exclusive_regions": regions,
    }
    return sets, summary


def compare_venn_membership(
    generated: pd.DataFrame,
    release: pd.DataFrame,
    observed_category: np.ndarray,
    expected_category: np.ndarray,
    contract: TableContract,
    generated_path: Path,
) -> dict[str, object] | None:
    if not contract.venn_conditions:
        return None
    observed_sets, observed_summary = venn_membership_summary(
        generated, observed_category, contract
    )
    expected_sets, expected_summary = venn_membership_summary(
        release, expected_category, contract
    )
    for condition in contract.venn_conditions:
        if observed_sets[condition] != expected_sets[condition]:
            difference = sorted(observed_sets[condition] ^ expected_sets[condition])
            raise ReleaseVerificationError(
                f"Venn membership differs for {generated_path}, condition={condition}, "
                f"direction={contract.venn_direction}: first differing gene={difference[0]}"
            )
    if observed_summary != expected_summary:
        raise ReleaseVerificationError(
            f"Venn-region membership differs for {generated_path}"
        )
    return observed_summary


def direction_label(value: float) -> str:
    if not np.isfinite(value):
        return "not_estimable"
    if value > 0:
        return "positive"
    if value < 0:
        return "negative"
    return "zero"


def compare_prespecified_interaction_genes(
    generated: pd.DataFrame,
    release: pd.DataFrame,
    parsed: dict[str, tuple[np.ndarray, np.ndarray]],
    contract: TableContract,
    generated_path: Path,
) -> list[dict[str, object]]:
    if not contract.prespecified_genes:
        return []
    observed_effect, expected_effect = parsed[contract.effect_column]
    observed_padj, expected_padj = parsed[contract.adjusted_p_column]
    records = []
    for gene in contract.prespecified_genes:
        rows = np.flatnonzero(generated["Gene_Symbol"].eq(gene).to_numpy())
        if len(rows) != 1:
            raise ReleaseVerificationError(
                f"Expected one {gene} row in interaction table {generated_path}"
            )
        row = int(rows[0])
        observed_conclusion = {
            "significant_at_fdr_0.05": bool(observed_padj[row] < 0.05),
            "direction": direction_label(observed_effect[row]),
        }
        expected_conclusion = {
            "significant_at_fdr_0.05": bool(expected_padj[row] < 0.05),
            "direction": direction_label(expected_effect[row]),
        }
        if observed_conclusion != expected_conclusion:
            raise ReleaseVerificationError(
                f"Prespecified interaction conclusion differs for {gene} in "
                f"{generated_path}: regenerated={observed_conclusion!r}, "
                f"release={expected_conclusion!r}"
            )
        records.append(
            {
                "gene": gene,
                **observed_conclusion,
                "generated_log2_fold_change": (
                    float(observed_effect[row])
                    if np.isfinite(observed_effect[row])
                    else None
                ),
                "release_log2_fold_change": (
                    float(expected_effect[row])
                    if np.isfinite(expected_effect[row])
                    else None
                ),
                "generated_adjusted_p": (
                    float(observed_padj[row])
                    if np.isfinite(observed_padj[row])
                    else None
                ),
                "release_adjusted_p": (
                    float(expected_padj[row])
                    if np.isfinite(expected_padj[row])
                    else None
                ),
            }
        )
    return records


def compare_guarded_figure_label_effects(
    generated: pd.DataFrame,
    release: pd.DataFrame,
    parsed: dict[str, tuple[np.ndarray, np.ndarray]],
    contract: TableContract,
    generated_path: Path,
) -> dict[str, object] | None:
    """Keep displayed label-gene fold changes stable to manuscript precision."""

    if not contract.guarded_effect_genes:
        return None
    if "condition" not in generated or "gene_symbol" not in generated:
        raise ReleaseVerificationError(
            f"Label-gene effect guard requires condition and gene_symbol in {generated_path}"
        )
    conditions = tuple(generated["condition"].drop_duplicates().tolist())
    observed_effect, expected_effect = parsed[contract.effect_column]
    records: list[dict[str, object]] = []
    maximum_delta = 0.0
    maximum_key: dict[str, str] | None = None
    estimable_rows = 0
    for condition in conditions:
        for gene in contract.guarded_effect_genes:
            selected = (
                generated["condition"].eq(condition)
                & generated["gene_symbol"].eq(gene)
            ).to_numpy()
            rows = np.flatnonzero(selected)
            if len(rows) != 1:
                raise ReleaseVerificationError(
                    f"Expected one label-gene row for condition={condition}, "
                    f"gene={gene} in {generated_path}"
                )
            row = int(rows[0])
            observed = observed_effect[row]
            expected = expected_effect[row]
            if np.isfinite(observed):
                estimable_rows += 1
                absolute_delta: float | None = float(abs(observed - expected))
                if maximum_key is None or absolute_delta > maximum_delta:
                    maximum_delta = absolute_delta
                    maximum_key = key_at(generated, row, contract)
                if absolute_delta > MANUSCRIPT_LABEL_LFC_ATOL:
                    raise ReleaseVerificationError(
                        f"Manuscript label-gene log2FC differs by more than "
                        f"{MANUSCRIPT_LABEL_LFC_ATOL:g} in {generated_path}, "
                        f"key={key_at(generated, row, contract)!r}: "
                        f"regenerated={observed:.17g}, release={expected:.17g}, "
                        f"absolute_delta={absolute_delta:.6g}"
                    )
                observed_value: float | None = float(observed)
                expected_value: float | None = float(expected)
            else:
                absolute_delta = None
                observed_value = None
                expected_value = None
            records.append(
                {
                    "key": key_at(generated, row, contract),
                    "generated_log2_fold_change": observed_value,
                    "release_log2_fold_change": expected_value,
                    "absolute_delta": absolute_delta,
                }
            )
    expected_rows = len(conditions) * len(contract.guarded_effect_genes)
    if len(records) != expected_rows:
        raise ReleaseVerificationError(
            f"Incomplete manuscript label-gene guard in {generated_path}"
        )
    return {
        "absolute_tolerance": MANUSCRIPT_LABEL_LFC_ATOL,
        "conditions": list(conditions),
        "genes": list(contract.guarded_effect_genes),
        "expected_rows": expected_rows,
        "estimable_rows": estimable_rows,
        "not_estimable_rows": expected_rows - estimable_rows,
        "max_absolute_delta": maximum_delta,
        "key_at_max_absolute_delta": maximum_key,
        "all_within_tolerance": True,
        "rows": records,
    }


def compare_scientific_outcomes(
    generated: pd.DataFrame,
    release: pd.DataFrame,
    parsed: dict[str, tuple[np.ndarray, np.ndarray]],
    contract: TableContract,
    generated_path: Path,
) -> dict[str, object]:
    observed_base, expected_base = parsed[contract.base_mean_column]
    observed_effect, expected_effect = parsed[contract.effect_column]
    observed_padj, expected_padj = parsed[contract.adjusted_p_column]

    base_masks = (observed_base >= 30.0, expected_base >= 30.0)
    up_masks = (observed_effect > 1.0, expected_effect > 1.0)
    down_masks = (observed_effect < -1.0, expected_effect < -1.0)
    require_equal_mask(
        *base_masks,
        "baseMean >= 30",
        generated,
        contract,
        observed_base,
        expected_base,
        generated_path,
    )
    require_equal_mask(
        *up_masks,
        "log2FC > 1",
        generated,
        contract,
        observed_effect,
        expected_effect,
        generated_path,
    )
    require_equal_mask(
        *down_masks,
        "log2FC < -1",
        generated,
        contract,
        observed_effect,
        expected_effect,
        generated_path,
    )

    observed_category = combined_category(
        observed_base, observed_effect, observed_padj
    )
    expected_category = combined_category(
        expected_base, expected_effect, expected_padj
    )
    require_equal_mask(
        observed_category,
        expected_category,
        "combined DEG category (baseMean >= 30, adjusted p < 0.05, |log2FC| > 1)",
        generated,
        contract,
        observed_effect,
        expected_effect,
        generated_path,
    )
    padj_flips = standalone_padj_flips(
        generated,
        contract,
        observed_padj,
        expected_padj,
        observed_category,
        expected_category,
    )
    if not padj_flips["outcome_categories_unchanged"]:
        raise ReleaseVerificationError(
            f"An adjusted-p threshold flip changes a combined DEG category in {generated_path}"
        )

    venn = compare_venn_membership(
        generated,
        release,
        observed_category,
        expected_category,
        contract,
        generated_path,
    )
    prespecified = compare_prespecified_interaction_genes(
        generated,
        release,
        parsed,
        contract,
        generated_path,
    )
    label_gene_effects = compare_guarded_figure_label_effects(
        generated,
        release,
        parsed,
        contract,
        generated_path,
    )
    return {
        "exact_threshold_masks": {
            "base_mean_ge_30": int(base_masks[0].sum()),
            "log2_fold_change_gt_1": int(up_masks[0].sum()),
            "log2_fold_change_lt_minus_1": int(down_masks[0].sum()),
        },
        "exact_combined_deg_categories": {
            "up": int((observed_category == 1).sum()),
            "down": int((observed_category == -1).sum()),
            "not_deg": int((observed_category == 0).sum()),
        },
        "standalone_adjusted_p_threshold_flips": padj_flips,
        "venn_membership": venn,
        "prespecified_interaction_genes": prespecified,
        "manuscript_label_gene_effect_guard": label_gene_effects,
    }


def compare_table(
    generated_path: Path,
    release_path: Path,
    contract: TableContract,
) -> dict[str, object]:
    generated = load_text_table(generated_path)
    release = load_text_table(release_path)
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
        validate_integer_columns(frame, contract.integer_columns, path)

    exact_columns = contract.key_columns + contract.exact_columns
    difference = first_exact_difference(generated, release, exact_columns)
    if difference is not None:
        row, column = difference
        raise ReleaseVerificationError(
            f"Exact semantic field differs for {generated_path}, "
            f"key={key_at(generated, row, contract)!r}, column={column}: "
            f"regenerated={generated.at[row, column]!r}, "
            f"release={release.at[row, column]!r}"
        )

    parsed, numeric_deltas = parse_and_compare_numeric_columns(
        generated, release, contract, generated_path, release_path
    )
    observed_values = {column: values[0] for column, values in parsed.items()}
    expected_values = {column: values[1] for column, values in parsed.items()}
    wald_identity = {
        "generated": validate_wald_identity(
            generated, observed_values, contract, generated_path
        ),
        "release": validate_wald_identity(
            release, expected_values, contract, release_path
        ),
    }
    outcomes = compare_scientific_outcomes(
        generated, release, parsed, contract, generated_path
    )

    generated_compressed_sha256 = sha256_file(generated_path)
    release_compressed_sha256 = sha256_file(release_path)
    generated_payload_sha256 = sha256_gzip_payload(generated_path)
    release_payload_sha256 = sha256_gzip_payload(release_path)
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
        "exact_schema_keys_order_text_integer_fields_and_na_masks": True,
        "wald_identity": wald_identity,
        "numeric_deltas_diagnostic_only": numeric_deltas,
        "scientific_outcomes": outcomes,
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
            "hashes_equal": generated.at[row, "sha256"] == release.at[row, "sha256"],
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


def new_report(generated_dir: Path, release_dir: Path) -> dict[str, object]:
    return {
        "schema_version": 2,
        "generated_dir": str(generated_dir),
        "release_dir": str(release_dir),
        "verification_policy": {
            "raw_numeric_deltas": "reported, not globally thresholded",
            "outcome_categories": "exact",
            "wald_identity_rtol": WALD_RTOL,
            "wald_identity_atol": WALD_ATOL,
            "manuscript_label_gene_lfc_atol": MANUSCRIPT_LABEL_LFC_ATOL,
        },
        "status": "failed",
        "exact_checks": {
            "manifest_excluding_output_sha256": "not_started",
            "generated_output_sha256_self_consistency": "not_started",
            "run_metadata": "not_started",
            "environment_freeze": "not_started",
        },
        "table_checks": {
            contract.release_filename: "not_started" for contract in TABLE_CONTRACTS
        },
        "tables": [],
    }


def run_stage(progress: dict[str, str], name: str, function, *args):
    progress[name] = "in_progress"
    try:
        result = function(*args)
    except Exception:
        progress[name] = "failed"
        raise
    progress[name] = "passed"
    return result


def verify_release(
    generated_dir: Path,
    release_dir: Path,
    report: dict[str, object] | None = None,
) -> list[dict[str, object]]:
    if report is None:
        report = new_report(generated_dir, release_dir)
    exact_checks = report["exact_checks"]
    table_checks = report["table_checks"]
    assert isinstance(exact_checks, dict)
    assert isinstance(table_checks, dict)
    manifest = run_stage(
        exact_checks,
        "manifest_excluding_output_sha256",
        compare_manifest,
        generated_dir / "analysis_manifest.tsv",
        release_dir / "analysis_manifest.tsv",
    )
    run_stage(
        exact_checks,
        "generated_output_sha256_self_consistency",
        verify_complete_exports,
        generated_dir,
        manifest,
    )
    run_stage(
        exact_checks,
        "run_metadata",
        compare_metadata,
        generated_dir / "run_metadata.json",
        release_dir / "run_metadata.json",
    )
    run_stage(
        exact_checks,
        "environment_freeze",
        compare_environment,
        generated_dir / "environment.freeze.txt",
        release_dir / "environment.freeze.txt",
    )

    summaries: list[dict[str, object]] = []
    for contract in TABLE_CONTRACTS:
        summary = run_stage(
            table_checks,
            contract.release_filename,
            compare_table,
            generated_dir / contract.generated_relative_path,
            release_dir / contract.release_filename,
            contract,
        )
        summaries.append(summary)
        report["tables"] = summaries.copy()
    report["manifest_output_hash_provenance"] = manifest_hash_provenance(
        generated_dir / "analysis_manifest.tsv",
        release_dir / "analysis_manifest.tsv",
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
    report = new_report(args.generated_dir, args.release_dir)
    error: Exception | None = None
    summaries: list[dict[str, object]] = []
    try:
        summaries = verify_release(args.generated_dir, args.release_dir, report)
        report["status"] = "passed"
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
        outcome = summary["scientific_outcomes"]
        categories = outcome["exact_combined_deg_categories"]
        flips = outcome["standalone_adjusted_p_threshold_flips"]["count"]
        print(
            "OUTCOME_MATCH "
            f"{summary['file']}: rows={summary['rows']}, "
            f"up={categories['up']}, down={categories['down']}, "
            f"standalone_padj_flips={flips}"
        )
    print(
        "BULK_RNASEQ_RESULTS_OK "
        "(exact structure, missingness, threshold outcomes, Venn membership, "
        "label-gene fold changes, prespecified interaction conclusions, and "
        "internal Wald identities)"
    )


if __name__ == "__main__":
    main()
