#!/usr/bin/env python3
"""Verify publication-facing bulk RNA-seq figure files and numeric contracts."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import struct
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
EXPECTED_DPI = 600.0
EXPECTED_PNG_PIXELS_PER_METRE = round(EXPECTED_DPI / 0.0254)
PUBLICATION_CONDITION_ORDER = ("TNF_IFNG", "IFNG", "TNF")
P_VALUE_FLOOR = 1e-300
DISPLAYED_Y_CAP = 300.0
WALD_95_Z = 1.959963984540054
FIGURE_2C_GENES = ("ICAM1", "IRF1")
FIGURE_2C_CONDITION_DISPLAY = {
    "TNF_IFNG": "TNF + IFNγ",
    "IFNG": "IFNγ",
    "TNF": "TNF",
}
FIGURE_2C_TSV_COLUMNS = (
    "condition",
    "condition_display",
    "gene_symbol",
    "contrast",
    "base_mean",
    "log2_fold_change_ko1_vs_wt",
    "lfc_se",
    "wald_ci_95_lower",
    "wald_ci_95_upper",
    "wald_statistic_ko1_vs_wt",
    "p_value",
    "adjusted_p_value_bh",
    "significance_code",
)


@dataclass(frozen=True)
class FigureSpec:
    name: str
    adapter: Path
    effect_column: str
    contract: Path
    figure_dir: Path
    triptych_stem: str
    individual_stems: tuple[str, ...]
    label_genes: tuple[str, ...]


def figure_specs(root: Path) -> tuple[FigureSpec, ...]:
    derived = root / "data" / "experimental" / "bulk_rnaseq" / "derived"
    return (
        FigureSpec(
            name="Figure 1B",
            adapter=derived
            / "figure_1b_1c_wt_cytokine_contrasts.unfiltered.tsv.gz",
            effect_column="log2_fold_change_treatment_vs_untreated",
            contract=root
            / "results"
            / "figure_1"
            / "Figure_1B_volcano_output_contract.tsv",
            figure_dir=root / "figures" / "figure_1",
            triptych_stem="Figure_1B_triptych",
            individual_stems=(
                "Figure_1B_TNF_IFNg_volcano",
                "Figure_1B_IFNg_volcano",
                "Figure_1B_TNF_volcano",
            ),
            label_genes=(
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
            ),
        ),
        FigureSpec(
            name="Figure 2B",
            adapter=derived
            / "figure_2b_s2d_tnfr1_ko1_vs_wt_matched_treatments.unfiltered.tsv.gz",
            effect_column="log2_fold_change_ko1_vs_wt",
            contract=root
            / "results"
            / "figure_2"
            / "Figure_2B_volcano_output_contract.tsv",
            figure_dir=root / "figures" / "figure_2",
            triptych_stem="Figure_2B_triptych",
            individual_stems=(
                "Figure_2B_TNF_IFNg",
                "Figure_2B_IFNg",
                "Figure_2B_TNF",
            ),
            label_genes=("ICAM1", "MLKL", "GSDME", "IRF1"),
        ),
    )


def parse_finite(value: str) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def figure_2c_significance_code(adjusted_p: float) -> str:
    if adjusted_p < 0.0001:
        return "****"
    if adjusted_p < 0.001:
        return "***"
    if adjusted_p < 0.01:
        return "**"
    if adjusted_p < 0.05:
        return "*"
    return "ns"


def expected_figure_2c_rows(root: Path) -> list[dict[str, object]]:
    adapter = (
        root
        / "data"
        / "experimental"
        / "bulk_rnaseq"
        / "derived"
        / "figure_2b_s2d_tnfr1_ko1_vs_wt_matched_treatments.unfiltered.tsv.gz"
    )
    selected: dict[tuple[str, str], dict[str, object]] = {}
    with gzip.open(adapter, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "condition",
            "gene_symbol",
            "base_mean",
            "log2_fold_change_ko1_vs_wt",
            "lfc_se",
            "wald_statistic_ko1_vs_wt",
            "p_value",
            "adjusted_p_value",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{adapter}: missing Figure 2C columns")
        for source in reader:
            condition = source["condition"]
            gene = source["gene_symbol"].strip().upper()
            if (
                condition not in PUBLICATION_CONDITION_ORDER
                or gene not in FIGURE_2C_GENES
            ):
                continue
            key = (condition, gene)
            if key in selected:
                raise ValueError(f"{adapter}: duplicate Figure 2C row {key}")
            numeric = {
                field: parse_finite(source[field])
                for field in required
                if field not in {"condition", "gene_symbol"}
            }
            if any(value is None for value in numeric.values()):
                raise ValueError(f"{adapter}: incomplete Figure 2C row {key}")
            selected[key] = numeric

    expected_keys = [
        (condition, gene)
        for condition in PUBLICATION_CONDITION_ORDER
        for gene in FIGURE_2C_GENES
    ]
    if set(selected) != set(expected_keys):
        raise ValueError(f"{adapter}: expected exactly six Figure 2C rows")

    rows: list[dict[str, object]] = []
    for condition, gene in expected_keys:
        source = selected[(condition, gene)]
        effect = float(source["log2_fold_change_ko1_vs_wt"])
        standard_error = float(source["lfc_se"])
        wald_statistic = float(source["wald_statistic_ko1_vs_wt"])
        adjusted_p = float(source["adjusted_p_value"])
        if standard_error <= 0 or not math.isclose(
            effect / standard_error,
            wald_statistic,
            rel_tol=1e-6,
            abs_tol=1e-8,
        ):
            raise ValueError(
                f"{adapter}: inconsistent Figure 2C Wald triplet for "
                f"{condition}/{gene}"
            )
        if not 0 <= adjusted_p <= 1:
            raise ValueError(
                f"{adapter}: invalid adjusted P value for {condition}/{gene}"
            )
        rows.append(
            {
                "condition": condition,
                "condition_display": FIGURE_2C_CONDITION_DISPLAY[condition],
                "gene_symbol": gene,
                "contrast": "TNFR1-KO1 vs WT within treatment",
                "base_mean": float(source["base_mean"]),
                "log2_fold_change_ko1_vs_wt": effect,
                "lfc_se": standard_error,
                "wald_ci_95_lower": effect - WALD_95_Z * standard_error,
                "wald_ci_95_upper": effect + WALD_95_Z * standard_error,
                "wald_statistic_ko1_vs_wt": wald_statistic,
                "p_value": float(source["p_value"]),
                "adjusted_p_value_bh": adjusted_p,
                "significance_code": figure_2c_significance_code(adjusted_p),
            }
        )
    return rows


def expected_contract_rows(spec: FigureSpec) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, object]]] = {
        condition: [] for condition in PUBLICATION_CONDITION_ORDER
    }
    with gzip.open(spec.adapter, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "condition",
            "gene_symbol",
            "base_mean",
            spec.effect_column,
            "adjusted_p_value",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{spec.adapter}: missing required columns")
        for source in reader:
            condition = source["condition"]
            if condition not in grouped:
                continue
            effect = parse_finite(source[spec.effect_column])
            adjusted_p = parse_finite(source["adjusted_p_value"])
            if effect is None or adjusted_p is None:
                continue
            base_mean = parse_finite(source["base_mean"])
            grouped[condition].append(
                {
                    "gene_symbol": source["gene_symbol"].strip().upper(),
                    "base_mean": base_mean,
                    "effect": effect,
                    "adjusted_p_value": adjusted_p,
                }
            )

    all_effects = [
        float(row["effect"])
        for condition in PUBLICATION_CONDITION_ORDER
        for row in grouped[condition]
    ]
    if not all_effects:
        raise ValueError(f"{spec.adapter}: no plottable effects")
    outward_limit = max(1, math.ceil(max(abs(value) for value in all_effects)))

    expected: list[dict[str, object]] = []
    for panel_order, condition in enumerate(PUBLICATION_CONDITION_ORDER, start=1):
        rows = grouped[condition]
        if not rows:
            raise ValueError(f"{spec.adapter}: no plottable rows for {condition}")
        effects = [float(row["effect"]) for row in rows]
        adjusted = [float(row["adjusted_p_value"]) for row in rows]
        finite_uncapped = [-math.log10(value) for value in adjusted if value > 0]
        observed_symbols = {
            str(row["gene_symbol"])
            for row in rows
            if row["base_mean"] is not None and float(row["base_mean"]) >= 30
        }
        plotted_labels = [gene for gene in spec.label_genes if gene in observed_symbols]
        expected.append(
            {
                "panel_order": panel_order,
                "condition": condition,
                "plottable_modelled_symbols": len(rows),
                "minimum_log2_fold_change": min(effects),
                "maximum_log2_fold_change": max(effects),
                "maximum_absolute_log2_fold_change": max(abs(value) for value in effects),
                "common_x_min": -outward_limit,
                "common_x_max": outward_limit,
                "adjusted_p_value_floor": P_VALUE_FLOOR,
                "displayed_y_cap": DISPLAYED_Y_CAP,
                "points_at_displayed_y_cap": sum(
                    value <= P_VALUE_FLOOR for value in adjusted
                ),
                "maximum_finite_uncapped_minus_log10_adjusted_p": max(
                    finite_uncapped
                ),
                "requested_label_count": len(spec.label_genes),
                "plotted_label_count": len(plotted_labels),
                "plotted_label_genes": ";".join(plotted_labels),
            }
        )
    return expected


def assert_close(
    observed: str,
    expected: float,
    context: str,
    *,
    abs_tol: float = 1e-12,
) -> None:
    value = float(observed)
    if not math.isclose(value, expected, rel_tol=1e-10, abs_tol=abs_tol):
        raise ValueError(f"{context}: observed {value!r}, expected {expected!r}")


def verify_contract(spec: FigureSpec) -> None:
    expected = expected_contract_rows(spec)
    with spec.contract.open("r", encoding="utf-8", newline="") as handle:
        observed = list(csv.DictReader(handle, delimiter="\t"))
    if len(observed) != len(expected):
        raise ValueError(
            f"{spec.contract}: observed {len(observed)} rows, expected {len(expected)}"
        )

    integer_fields = (
        "panel_order",
        "plottable_modelled_symbols",
        "points_at_displayed_y_cap",
        "requested_label_count",
        "plotted_label_count",
    )
    numeric_fields = (
        "minimum_log2_fold_change",
        "maximum_log2_fold_change",
        "maximum_absolute_log2_fold_change",
        "common_x_min",
        "common_x_max",
        "adjusted_p_value_floor",
        "displayed_y_cap",
        "maximum_finite_uncapped_minus_log10_adjusted_p",
    )
    text_fields = ("condition", "plotted_label_genes")
    for index, (actual, wanted) in enumerate(zip(observed, expected), start=1):
        context = f"{spec.contract}, row {index}"
        for field in integer_fields:
            if int(actual[field]) != wanted[field]:
                raise ValueError(
                    f"{context}, {field}: observed {actual[field]!r}, "
                    f"expected {wanted[field]!r}"
                )
        for field in numeric_fields:
            assert_close(
                actual[field],
                float(wanted[field]),
                f"{context}, {field}",
                abs_tol=0.0 if field == "adjusted_p_value_floor" else 1e-12,
            )
        for field in text_fields:
            if actual[field] != wanted[field]:
                raise ValueError(
                    f"{context}, {field}: observed {actual[field]!r}, "
                    f"expected {wanted[field]!r}"
                )


def read_png_metadata(path: Path) -> dict[str, object]:
    data = path.read_bytes()
    if not data.startswith(PNG_SIGNATURE):
        raise ValueError(f"{path}: not a PNG file")
    offset = len(PNG_SIGNATURE)
    width = height = None
    pixels_per_metre = None
    while offset + 12 <= len(data):
        length = struct.unpack(">I", data[offset : offset + 4])[0]
        chunk_type = data[offset + 4 : offset + 8]
        payload = data[offset + 8 : offset + 8 + length]
        if len(payload) != length:
            raise ValueError(f"{path}: truncated PNG chunk")
        if chunk_type == b"IHDR":
            width, height = struct.unpack(">II", payload[:8])
        elif chunk_type == b"pHYs" and length == 9:
            x_ppm, y_ppm, unit = struct.unpack(">IIB", payload)
            pixels_per_metre = (x_ppm, y_ppm, unit)
        offset += 12 + length
        if chunk_type == b"IEND":
            break
    if width is None or height is None:
        raise ValueError(f"{path}: PNG IHDR is missing")
    return {
        "width": width,
        "height": height,
        "pixels_per_metre": pixels_per_metre,
    }


TIFF_TYPE_SIZES = {1: 1, 2: 1, 3: 2, 4: 4, 5: 8}


def read_tiff_tags(path: Path) -> dict[int, tuple[float, ...]]:
    data = path.read_bytes()
    if len(data) < 8 or data[:2] not in (b"II", b"MM"):
        raise ValueError(f"{path}: not a classic TIFF file")
    endian = "<" if data[:2] == b"II" else ">"
    if struct.unpack(endian + "H", data[2:4])[0] != 42:
        raise ValueError(f"{path}: unsupported TIFF variant")
    ifd_offset = struct.unpack(endian + "I", data[4:8])[0]
    if ifd_offset + 2 > len(data):
        raise ValueError(f"{path}: invalid TIFF IFD offset")
    count = struct.unpack(endian + "H", data[ifd_offset : ifd_offset + 2])[0]
    tags: dict[int, tuple[float, ...]] = {}
    entry_offset = ifd_offset + 2
    for _ in range(count):
        entry = data[entry_offset : entry_offset + 12]
        if len(entry) != 12:
            raise ValueError(f"{path}: truncated TIFF IFD")
        tag, value_type, value_count = struct.unpack(endian + "HHI", entry[:8])
        size = TIFF_TYPE_SIZES.get(value_type)
        if size is None:
            entry_offset += 12
            continue
        total_size = size * value_count
        if total_size <= 4:
            raw = entry[8 : 8 + total_size]
        else:
            value_offset = struct.unpack(endian + "I", entry[8:12])[0]
            raw = data[value_offset : value_offset + total_size]
        if len(raw) != total_size:
            raise ValueError(f"{path}: truncated TIFF tag {tag}")
        if value_type == 3:
            values = struct.unpack(endian + "H" * value_count, raw)
            tags[tag] = tuple(float(value) for value in values)
        elif value_type == 4:
            values = struct.unpack(endian + "I" * value_count, raw)
            tags[tag] = tuple(float(value) for value in values)
        elif value_type == 5:
            integers = struct.unpack(endian + "I" * (2 * value_count), raw)
            tags[tag] = tuple(
                integers[index] / integers[index + 1]
                for index in range(0, len(integers), 2)
            )
        entry_offset += 12
    return tags


def verify_png(path: Path, expected_pixels: tuple[int, int]) -> None:
    metadata = read_png_metadata(path)
    observed_pixels = (metadata["width"], metadata["height"])
    if observed_pixels != expected_pixels:
        raise ValueError(
            f"{path}: observed {observed_pixels} pixels, expected {expected_pixels}"
        )
    density = metadata["pixels_per_metre"]
    if density is None:
        raise ValueError(f"{path}: PNG physical-resolution metadata is missing")
    x_ppm, y_ppm, unit = density
    if unit != 1 or abs(x_ppm - EXPECTED_PNG_PIXELS_PER_METRE) > 1 or abs(
        y_ppm - EXPECTED_PNG_PIXELS_PER_METRE
    ) > 1:
        raise ValueError(f"{path}: PNG physical resolution is not 600 dpi")


def scalar_tag(tags: dict[int, tuple[float, ...]], tag: int, path: Path) -> float:
    values = tags.get(tag)
    if values is None or len(values) != 1:
        raise ValueError(f"{path}: required scalar TIFF tag {tag} is missing")
    return values[0]


def verify_tiff(path: Path, expected_pixels: tuple[int, int]) -> None:
    tags = read_tiff_tags(path)
    observed_pixels = (
        int(scalar_tag(tags, 256, path)),
        int(scalar_tag(tags, 257, path)),
    )
    if observed_pixels != expected_pixels:
        raise ValueError(
            f"{path}: observed {observed_pixels} pixels, expected {expected_pixels}"
        )
    if int(scalar_tag(tags, 259, path)) != 5:
        raise ValueError(f"{path}: TIFF compression is not LZW")
    if int(scalar_tag(tags, 296, path)) != 2:
        raise ValueError(f"{path}: TIFF resolution unit is not inches")
    for tag in (282, 283):
        if not math.isclose(
            scalar_tag(tags, tag, path), EXPECTED_DPI, rel_tol=0, abs_tol=0.1
        ):
            raise ValueError(f"{path}: TIFF resolution is not 600 dpi")


def verify_image_pair(stem: Path, expected_pixels: tuple[int, int]) -> None:
    verify_png(stem.with_suffix(".png"), expected_pixels)
    verify_tiff(stem.with_suffix(".tiff"), expected_pixels)


def verify_figure(spec: FigureSpec) -> None:
    verify_contract(spec)
    verify_image_pair(spec.figure_dir / spec.triptych_stem, (3900, 1800))
    for stem in spec.individual_stems:
        verify_image_pair(spec.figure_dir / stem, (4800, 3600))


def verify_figure_2c(root: Path) -> None:
    table_path = (
        root / "results" / "figure_2" / "Figure_2C_ICAM1_IRF1_effects.tsv"
    )
    with table_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != FIGURE_2C_TSV_COLUMNS:
            raise ValueError(f"{table_path}: unexpected columns")
        observed = list(reader)
    expected = expected_figure_2c_rows(root)
    if len(observed) != len(expected):
        raise ValueError(
            f"{table_path}: observed {len(observed)} rows, expected {len(expected)}"
        )

    text_fields = (
        "condition",
        "condition_display",
        "gene_symbol",
        "contrast",
        "significance_code",
    )
    numeric_fields = tuple(
        field for field in FIGURE_2C_TSV_COLUMNS if field not in text_fields
    )
    for index, (actual, wanted) in enumerate(zip(observed, expected), start=1):
        context = f"{table_path}, row {index}"
        for field in text_fields:
            if actual[field] != wanted[field]:
                raise ValueError(
                    f"{context}, {field}: observed {actual[field]!r}, "
                    f"expected {wanted[field]!r}"
                )
        for field in numeric_fields:
            assert_close(
                actual[field],
                float(wanted[field]),
                f"{context}, {field}",
                abs_tol=0.0
                if field in {"p_value", "adjusted_p_value_bh"}
                else 1e-12,
            )

    verify_image_pair(
        root / "figures" / "figure_2" / "Figure_2C_ICAM1_IRF1_effects",
        (4500, 2280),
    )


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repository-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    args = parser.parse_args(argv)
    root = args.repository_root.resolve()
    specs = figure_specs(root)
    for spec in specs:
        verify_figure(spec)
        expected = expected_contract_rows(spec)
        common_x = (expected[0]["common_x_min"], expected[0]["common_x_max"])
        capped = [row["points_at_displayed_y_cap"] for row in expected]
        print(
            f"{spec.name}: passed; common x range {common_x}; "
            f"points at y=300 by panel {capped}"
        )
    verify_figure_2c(root)
    print("Figure 2C: passed; six PyDESeq2 effects with 95% Wald CIs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
