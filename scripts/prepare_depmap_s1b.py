#!/usr/bin/env python3
"""Prepare the checksum-pinned DepMap 25Q2 table for Supplementary Figure S1B."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any


EXPRESSION_FILE_NAME = "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
MODEL_FILE_NAME = "Model.csv"

EXPECTED_EXPRESSION_SIZE = 521_526_607
EXPECTED_EXPRESSION_SHA256 = (
    "e0326e16eb23bea1be980fce315acb36b224dedd7af6b47e0ba37e7747dbcc47"
)
EXPECTED_MODEL_SIZE = 694_278
EXPECTED_MODEL_SHA256 = (
    "b096e03bfefdc2679211545ddbf1bb7878d69ffde07ae335af5b968a7883733c"
)

OUTPUT_COLUMNS = (
    "ProfileID",
    "ModelID",
    "OncotreePrimaryDisease",
    "RIPK3_log2_TPM_plus_1",
    "NLRP3_log2_TPM_plus_1",
    "RIPK3_below_threshold",
    "NLRP3_below_threshold",
    "threshold_category",
)


def file_digest(path: Path, algorithm: str = "sha256") -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stream_digest(handle: Any, algorithm: str = "sha256") -> str:
    digest = hashlib.new(algorithm)
    for block in iter(lambda: handle.read(1024 * 1024), b""):
        digest.update(block)
    return digest.hexdigest()


def assert_source_file(
    path: Path,
    *,
    expected_size: int,
    expected_sha256: str,
    expected_md5: str | None = None,
) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing input: {path}")
    observed_size = path.stat().st_size
    if observed_size != expected_size:
        raise ValueError(
            f"Byte-size mismatch for {path}: expected {expected_size}, "
            f"observed {observed_size}"
        )
    observed_sha256 = file_digest(path)
    if observed_sha256.lower() != expected_sha256.lower():
        raise ValueError(
            f"SHA-256 mismatch for {path}: observed {observed_sha256}"
        )
    result: dict[str, Any] = {
        "filename": path.name,
        "size_bytes": observed_size,
        "sha256": observed_sha256,
    }
    if expected_md5 is not None:
        observed_md5 = file_digest(path, "md5")
        if observed_md5.lower() != expected_md5.lower():
            raise ValueError(f"MD5 mismatch for {path}: observed {observed_md5}")
        result["md5"] = observed_md5
    return result


def load_models(path: Path) -> tuple[dict[str, tuple[str, str]], dict[str, Any]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("Model.csv has no header")
        required = {"ModelID", "ModelType", "OncotreePrimaryDisease"}
        missing = sorted(required.difference(reader.fieldnames))
        if missing:
            raise ValueError(f"Model.csv is missing columns: {', '.join(missing)}")
        tissue_origin_present = "TissueOrigin" in reader.fieldnames
        models: dict[str, tuple[str, str]] = {}
        model_types: Counter[str] = Counter()
        non_cancerous = 0
        tissue_origin_nonempty = 0
        for row in reader:
            model_id = row["ModelID"].strip()
            model_type = row["ModelType"].strip()
            disease = row["OncotreePrimaryDisease"].strip()
            if not model_id or model_id in models:
                raise ValueError("Model.csv ModelID values must be non-empty and unique")
            if not disease:
                raise ValueError("OncotreePrimaryDisease must be non-empty")
            models[model_id] = (model_type, disease)
            model_types[model_type] += 1
            non_cancerous += disease.upper() == "NON-CANCEROUS"
            if tissue_origin_present:
                tissue_origin_nonempty += bool(row["TissueOrigin"].strip())

    qc = {
        "rows": len(models),
        "unique_model_ids": len(models),
        "model_type_counts": dict(sorted(model_types.items())),
        "non_cancerous_rows": non_cancerous,
        "tissue_origin_present": tissue_origin_present,
        "tissue_origin_nonempty_rows": (
            tissue_origin_nonempty if tissue_origin_present else None
        ),
        "tissue_origin_filter_applied": False,
    }
    return models, qc


def load_default_expression(
    expression_path: Path,
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    gene_columns = ("RIPK3 (11035)", "NLRP3 (114548)")
    with expression_path.open("r", encoding="utf-8-sig", newline="") as text_handle:
        reader = csv.reader(text_handle)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError("Expression matrix is empty") from exc
        missing = [column for column in gene_columns if column not in header]
        if missing:
            raise ValueError(
                f"Expression matrix is missing gene columns: {', '.join(missing)}"
            )

        # The direct 25Q2 matrix stores ACH ModelID values in its first column
        # and has no ProfileID/is_default_entry metadata. Older matrices used
        # explicit metadata columns, which remain supported here.
        model_id_column = "ModelID" if "ModelID" in header else header[0]
        profile_id_column = "ProfileID" if "ProfileID" in header else model_id_column
        default_flag_column = (
            "is_default_entry" if "is_default_entry" in header else None
        )
        index = {column: header.index(column) for column in gene_columns}
        model_id_index = header.index(model_id_column)
        profile_id_index = header.index(profile_id_column)
        default_flag_index = (
            header.index(default_flag_column) if default_flag_column else None
        )
        default_rows: list[dict[str, str]] = []
        profile_ids: set[str] = set()
        default_model_ids: set[str] = set()
        source_rows = 0
        nondefault_rows = 0
        for values in reader:
            if not values or not any(value.strip() for value in values):
                continue
            source_rows += 1
            if default_flag_index is not None:
                flag = values[default_flag_index].strip().lower()
                if flag not in {"true", "false"}:
                    raise ValueError("is_default_entry contains a non-binary value")
                if flag == "false":
                    nondefault_rows += 1
                    continue
            profile_id = values[profile_id_index].strip()
            model_id = values[model_id_index].strip()
            ripk3 = values[index["RIPK3 (11035)"]].strip()
            nlrp3 = values[index["NLRP3 (114548)"]].strip()
            if not profile_id or profile_id in profile_ids:
                raise ValueError("Profile identifiers must be non-empty and unique")
            if not model_id or model_id in default_model_ids:
                raise ValueError("ModelID values must be non-empty and unique")
            if not model_id.startswith("ACH-") or len(model_id) != 10:
                raise ValueError(f"Unexpected ModelID format: {model_id}")
            if not math.isfinite(float(ripk3)) or not math.isfinite(float(nlrp3)):
                raise ValueError("RIPK3 and NLRP3 must be finite in retained rows")
            profile_ids.add(profile_id)
            default_model_ids.add(model_id)
            default_rows.append(
                {
                    "ProfileID": profile_id,
                    "ModelID": model_id,
                    "RIPK3_log2_TPM_plus_1": ripk3,
                    "NLRP3_log2_TPM_plus_1": nlrp3,
                }
            )

    if source_rows == 0 or len(default_rows) == 0:
        raise ValueError("Expression matrix contains no usable default rows")

    expression_qc = {
        "source_profile_rows": source_rows,
        "default_profile_rows": len(default_rows),
        "nondefault_profile_rows": nondefault_rows,
        "unique_default_model_ids": len(default_model_ids),
        "unique_default_profile_ids": len(profile_ids),
        "model_id_source_column": model_id_column or "first unnamed column",
        "profile_id_source_column": profile_id_column or "first unnamed column",
        "default_entry_filter_applied": default_flag_column is not None,
    }
    return default_rows, expression_qc


def join_annotations(
    default_rows: list[dict[str, str]],
    models: dict[str, tuple[str, str]],
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    missing = sorted({row["ModelID"] for row in default_rows}.difference(models))
    if missing:
        raise ValueError(
            f"Model.csv is missing {len(missing)} default ModelID values; "
            f"first: {missing[0]}"
        )
    eligible: list[dict[str, str]] = []
    default_cell_lines = 0
    default_non_cancerous = 0
    for row in default_rows:
        model_type, disease = models[row["ModelID"]]
        default_cell_lines += model_type.upper() == "CELL LINE"
        default_non_cancerous += (
            model_type.upper() == "CELL LINE" and disease.upper() == "NON-CANCEROUS"
        )
        if model_type.upper() != "CELL LINE" or disease.upper() == "NON-CANCEROUS":
            continue
        ripk3_low = float(row["RIPK3_log2_TPM_plus_1"]) < 0.5
        nlrp3_low = float(row["NLRP3_log2_TPM_plus_1"]) < 0.5
        if ripk3_low and nlrp3_low:
            category = "Both below threshold"
        elif ripk3_low:
            category = "RIPK3 below threshold only"
        elif nlrp3_low:
            category = "NLRP3 below threshold only"
        else:
            category = "Neither below threshold"
        eligible.append(
            {
                "ProfileID": row["ProfileID"],
                "ModelID": row["ModelID"],
                "OncotreePrimaryDisease": disease,
                "RIPK3_log2_TPM_plus_1": row["RIPK3_log2_TPM_plus_1"],
                "NLRP3_log2_TPM_plus_1": row["NLRP3_log2_TPM_plus_1"],
                "RIPK3_below_threshold": str(ripk3_low).upper(),
                "NLRP3_below_threshold": str(nlrp3_low).upper(),
                "threshold_category": category,
            }
        )
    eligible.sort(key=lambda row: row["ModelID"])
    if not eligible:
        raise ValueError("No eligible cancer cell-line models remained after filtering")
    threshold_counts = Counter()
    for row in eligible:
        ripk3_low = float(row["RIPK3_log2_TPM_plus_1"]) < 0.5
        nlrp3_low = float(row["NLRP3_log2_TPM_plus_1"]) < 0.5
        threshold_counts["RIPK3_below_threshold"] += ripk3_low
        threshold_counts["NLRP3_below_threshold"] += nlrp3_low
        threshold_counts["both_below_threshold"] += ripk3_low and nlrp3_low
        threshold_counts["RIPK3_below_threshold_only"] += (
            ripk3_low and not nlrp3_low
        )
        threshold_counts["NLRP3_below_threshold_only"] += (
            not ripk3_low and nlrp3_low
        )
        threshold_counts["neither_below_threshold"] += (
            not ripk3_low and not nlrp3_low
        )
    join_qc = {
        "default_profiles_matched_to_metadata": len(default_rows),
        "default_profiles_annotated_cell_line": default_cell_lines,
        "default_cell_lines_annotated_non_cancerous": default_non_cancerous,
        "eligible_models": len(eligible),
        "threshold_log2_tpm_plus_1": 0.5,
        "threshold_counts": dict(threshold_counts),
    }
    return eligible, join_qc


def write_deterministic_tsv_gz(
    rows: list[dict[str, str]], output_path: Path
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output_path.name}.", dir=output_path.parent
    )
    os.close(descriptor)
    temporary_path = Path(temporary_name)
    try:
        with temporary_path.open("wb") as binary_handle:
            with gzip.GzipFile(
                filename="",
                mode="wb",
                compresslevel=9,
                fileobj=binary_handle,
                mtime=0,
            ) as gzip_handle:
                with io.TextIOWrapper(
                    gzip_handle, encoding="utf-8", newline=""
                ) as text_handle:
                    writer = csv.DictWriter(
                        text_handle,
                        fieldnames=OUTPUT_COLUMNS,
                        delimiter="\t",
                        lineterminator="\n",
                        extrasaction="raise",
                    )
                    writer.writeheader()
                    writer.writerows(rows)
        temporary_path.replace(output_path)
    finally:
        temporary_path.unlink(missing_ok=True)


def write_json_atomic(payload: dict[str, Any], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output_path.name}.", dir=output_path.parent
    )
    os.close(descriptor)
    temporary_path = Path(temporary_name)
    try:
        temporary_path.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        temporary_path.replace(output_path)
    finally:
        temporary_path.unlink(missing_ok=True)


def write_summary_csv(
    output_path: Path, threshold_counts: dict[str, int], denominator: int
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output_path.name}.", dir=output_path.parent
    )
    os.close(descriptor)
    temporary_path = Path(temporary_name)
    labels = {
        "RIPK3_below_threshold": "RIPK3 below threshold",
        "NLRP3_below_threshold": "NLRP3 below threshold",
        "both_below_threshold": "Both below threshold",
        "RIPK3_below_threshold_only": "RIPK3 below threshold only",
        "NLRP3_below_threshold_only": "NLRP3 below threshold only",
        "neither_below_threshold": "Neither below threshold",
    }
    try:
        with temporary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=(
                    "metric",
                    "n",
                    "denominator",
                    "percent",
                    "cutoff_log2_tpm_plus_1",
                ),
                lineterminator="\n",
            )
            writer.writeheader()
            for key in labels:
                count = threshold_counts.get(key, 0)
                writer.writerow(
                    {
                        "metric": labels[key],
                        "n": count,
                        "denominator": denominator,
                        "percent": f"{100 * count / denominator:.1f}",
                        "cutoff_log2_tpm_plus_1": "0.5",
                    }
                )
        temporary_path.replace(output_path)
    finally:
        temporary_path.unlink(missing_ok=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--expression-file",
        "--expression-archive",
        dest="expression_file",
        type=Path,
        default=Path("data/depmap/raw") / EXPRESSION_FILE_NAME,
        help="Direct DepMap 25Q2 expression CSV (old --expression-archive alias retained)",
    )
    parser.add_argument(
        "--model-file",
        type=Path,
        default=Path("data/depmap/raw") / MODEL_FILE_NAME,
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/analysis/depmap_s1b_eligible_models.tsv.gz"),
    )
    parser.add_argument(
        "--qc-output",
        type=Path,
        default=Path("data/analysis/depmap_s1b_preparation_qc.json"),
    )
    parser.add_argument(
        "--provenance-output",
        type=Path,
        default=Path("data/analysis/depmap_s1b_source_provenance.json"),
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=Path("reference_results/depmap_s1b_statistics.csv"),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    expression_source = assert_source_file(
        args.expression_file,
        expected_size=EXPECTED_EXPRESSION_SIZE,
        expected_sha256=EXPECTED_EXPRESSION_SHA256,
    )
    expression_source["filename"] = EXPRESSION_FILE_NAME
    model_source = assert_source_file(
        args.model_file,
        expected_size=EXPECTED_MODEL_SIZE,
        expected_sha256=EXPECTED_MODEL_SHA256,
    )
    model_source["filename"] = MODEL_FILE_NAME
    models, model_qc = load_models(args.model_file)
    default_rows, expression_qc = load_default_expression(args.expression_file)
    eligible, join_qc = join_annotations(default_rows, models)
    write_deterministic_tsv_gz(eligible, args.output)
    output_source = {
        "filename": args.output.name,
        "size_bytes": args.output.stat().st_size,
        "sha256": file_digest(args.output),
        "rows": len(eligible),
        "columns": list(OUTPUT_COLUMNS),
        "sort_order": ["ModelID"],
        "compression": "gzip level 9; mtime 0; no original filename",
    }
    qc = {
        "contract_version": 1,
        "expression_qc": expression_qc,
        "model_qc": model_qc,
        "join_qc": join_qc,
        "derived_file": output_source,
    }
    provenance = {
        "contract_version": 1,
        "expression_release": "DepMap Public 25Q2",
        "expression_release_identity_status": "confirmed",
        "model_release": None,
        "model_release_identity_status": "unverified",
        "release_pair_status": "unverified",
        "same_release_pair": None,
        "release_identity_note": (
            "The expression matrix is DepMap Public 25Q2. Model.csv is "
            "separately checksum-pinned, but its release identity must be "
            "confirmed before a same-release claim is made."
        ),
        "source_files": {
            "expression_csv": expression_source,
            "model_metadata": model_source,
        },
        "derived_file": output_source,
        "population_rule": (
            "is_default_entry == True; ModelType == Cell Line; non-empty "
            "OncotreePrimaryDisease != Non-Cancerous; finite RIPK3 and NLRP3"
        ),
        "tissue_origin_filter_applied": False,
    }
    write_json_atomic(qc, args.qc_output)
    write_json_atomic(provenance, args.provenance_output)
    write_summary_csv(
        args.summary_output,
        join_qc["threshold_counts"],
        join_qc["eligible_models"],
    )
    print(
        json.dumps(
            {
                "derived_file": output_source,
                "qc": str(args.qc_output),
                "provenance": str(args.provenance_output),
                "summary": str(args.summary_output),
            }
        )
    )


if __name__ == "__main__":
    main()
