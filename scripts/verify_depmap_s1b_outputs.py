#!/usr/bin/env python3
"""Verify the publication-facing DepMap Supplementary Figure S1B outputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable

from verify_bulk_figure_outputs import verify_image_pair


EXPECTED_PROVENANCE = {
    "release_pair_status": "confirmed",
    "expression_release": "DepMap Public 25Q2",
    "model_release": "DepMap Public 25Q2",
    "model_release_identity_status": "confirmed",
    "same_release_pair": "TRUE",
    "derived_file": "depmap_s1b_eligible_models.tsv.gz",
    "population": (
        "DepMap cell-line models with a non-missing OncoTree primary-disease "
        "label other than Non-Cancerous"
    ),
    "n_models": "1591",
    "cutoff_log2_tpm_plus_1": "0.5",
    "RIPK3_below_threshold": "1003",
    "NLRP3_below_threshold": "1172",
    "both_below_threshold": "749",
    "RIPK3_below_threshold_only": "254",
    "NLRP3_below_threshold_only": "423",
    "neither_below_threshold": "165",
    "tissue_origin_filter_applied": "FALSE",
}


def read_provenance(path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or set(rows[0]) != {"field", "value"}:
        raise ValueError(f"{path}: unexpected provenance schema")
    if any(set(row) != {"field", "value"} for row in rows):
        raise ValueError(f"{path}: inconsistent provenance rows")
    fields = [row["field"] for row in rows]
    if len(fields) != len(set(fields)):
        raise ValueError(f"{path}: duplicate provenance fields")
    return {row["field"]: row["value"] for row in rows}


def verify_outputs(root: Path) -> None:
    output_dir = root / "results" / "supplementary_S1B"
    stem = output_dir / "Supplementary_Figure_S1B"
    verify_image_pair(stem, (4260, 3840))

    provenance_path = output_dir / "Supplementary_Figure_S1B_runtime_provenance.tsv"
    observed = read_provenance(provenance_path)
    if observed != EXPECTED_PROVENANCE:
        missing = sorted(set(EXPECTED_PROVENANCE) - set(observed))
        extra = sorted(set(observed) - set(EXPECTED_PROVENANCE))
        changed = sorted(
            key
            for key in set(observed) & set(EXPECTED_PROVENANCE)
            if observed[key] != EXPECTED_PROVENANCE[key]
        )
        raise ValueError(
            f"{provenance_path}: provenance mismatch; "
            f"missing={missing}, extra={extra}, changed={changed}"
        )

    session_path = output_dir / "sessionInfo.txt"
    session = session_path.read_text(encoding="utf-8")
    required_session_tokens = ("R version 4.4.3", "ggplot2_", "data.table_")
    missing_tokens = [token for token in required_session_tokens if token not in session]
    if missing_tokens:
        raise ValueError(f"{session_path}: missing tokens {missing_tokens}")


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repository-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    args = parser.parse_args(argv)
    verify_outputs(args.repository_root.resolve())
    print(
        "Supplementary Figure S1B: passed; n=1,591; "
        "PNG/TIFF 600 dpi; TIFF LZW; provenance locked as a confirmed "
        "DepMap Public 25Q2 source pair"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
