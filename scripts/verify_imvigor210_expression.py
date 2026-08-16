#!/usr/bin/env python3
"""Stream-verify the IMvigor210 expression export without exposing its data.

The required compatibility digest is independent of CSV quoting, line endings
and harmless decimal text differences.  Version 1 is framed exactly as follows
(all text is UTF-8 and every displayed newline is one LF byte):

    IMVIGOR210_EXPRESSION_SEMANTIC_V1\n
    scale=<decimal scale>\n
    rows=<row count>\n
    columns=<column count>\n
    column_ids\n
    <UTF-8 byte length>:<column identifier>\n       (once per column)
    rows\n
    <UTF-8 byte length>:<row identifier>\t<i1>\t...\t<iN>\n

Each ``i`` is the signed base-10 integer obtained by parsing the CSV field as
``Decimal``, rounding with ``ROUND_HALF_UP`` to the requested fixed scale, and
multiplying by ``10**scale``.  Zero is always encoded as ``0``.  Ordered row
and column identifier digests use the magic line
``IMVIGOR210_ORDERED_IDS_V1\n`` followed by the same length-prefixed identifier
records.  This is the executable reference implementation of that framing.

Only aggregate dimensions and cryptographic digests are emitted.  Exceptions
and JSON reports never contain a sample identifier, feature identifier or
expression value.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import tempfile
from decimal import Decimal, DecimalException, ROUND_HALF_UP, localcontext
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "resources" / "IMvigor210_expression_semantic_contract_v1.json"
EXPRESSION_MAGIC = b"IMVIGOR210_EXPRESSION_SEMANTIC_V1\n"
ORDERED_IDS_MAGIC = b"IMVIGOR210_ORDERED_IDS_V1\n"


class VerificationError(RuntimeError):
    """Raised for a non-identifying contract or input verification failure."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, help="IMvigor210 expression CSV")
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument(
        "--report",
        type=Path,
        help="Optional path for a redacted JSON verification report",
    )
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="Run the embedded non-sensitive framing test and exit",
    )
    arguments = parser.parse_args()
    if not arguments.self_test and arguments.input is None:
        parser.error("--input is required unless --self-test is used")
    return arguments


def _sha256() -> Any:
    return hashlib.sha256()


def _identifier_record(identifier: str) -> bytes:
    encoded = identifier.encode("utf-8")
    return str(len(encoded)).encode("ascii") + b":" + encoded + b"\n"


def _row_prefix(identifier: str) -> bytes:
    encoded = identifier.encode("utf-8")
    return str(len(encoded)).encode("ascii") + b":" + encoded


def _contract_int(mapping: dict[str, Any], key: str) -> int:
    value = mapping.get(key)
    if isinstance(value, bool) or not isinstance(value, int):
        raise VerificationError(f"Contract field {key!r} must be an integer")
    return value


def _contract_digest(value: object, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64:
        raise VerificationError(f"Contract field {label!r} is not a SHA-256 digest")
    try:
        int(value, 16)
    except ValueError as exc:
        raise VerificationError(
            f"Contract field {label!r} is not a SHA-256 digest"
        ) from exc
    return value.lower()


def load_contract(path: Path) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as handle:
            contract = json.load(handle)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise VerificationError("Could not read the semantic contract") from exc
    if not isinstance(contract, dict):
        raise VerificationError("Semantic contract root must be an object")
    if contract.get("semantic_contract_version") != 1:
        raise VerificationError("Unsupported semantic contract version")

    shape = contract.get("shape")
    identifiers = contract.get("ordered_identifier_digests")
    semantics = contract.get("semantic_digests")
    numeric = contract.get("numeric_constraints")
    analysis = contract.get("analysis_canonicalization")
    if not all(
        isinstance(item, dict)
        for item in (shape, identifiers, semantics, numeric, analysis)
    ):
        raise VerificationError("Semantic contract is missing a required object")
    rows = _contract_int(shape, "rows")
    columns = _contract_int(shape, "columns")
    values = _contract_int(shape, "values")
    if rows <= 0 or columns <= 0 or values != rows * columns:
        raise VerificationError("Semantic contract shape is inconsistent")
    if _contract_int(numeric, "zero_count") < 0:
        raise VerificationError("Semantic contract zero count is invalid")
    _contract_digest(identifiers.get("row_ids_sha256"), "row_ids_sha256")
    _contract_digest(identifiers.get("column_ids_sha256"), "column_ids_sha256")

    required_scale = _contract_int(semantics, "required_scale")
    diagnostic_scales = semantics.get("diagnostic_scales")
    reference = semantics.get("reference_sha256")
    if not isinstance(diagnostic_scales, list) or not isinstance(reference, dict):
        raise VerificationError("Semantic digest scale configuration is invalid")
    scales = [required_scale, *diagnostic_scales]
    if (
        any(isinstance(scale, bool) or not isinstance(scale, int) for scale in scales)
        or any(scale < 0 or scale > 18 for scale in scales)
        or len(scales) != len(set(scales))
    ):
        raise VerificationError("Semantic digest scales are invalid")
    for scale in scales:
        _contract_digest(reference.get(f"fixed{scale}"), f"fixed{scale}")
    if (
        analysis.get("required") is not True
        or analysis.get("scale") != required_scale
        or analysis.get("rounding") != "ROUND_HALF_UP"
    ):
        raise VerificationError(
            "Semantic and analysis canonicalization policies are inconsistent"
        )
    return contract


def _initial_semantic_digest(scale: int, rows: int, columns: int) -> Any:
    digest = _sha256()
    digest.update(EXPRESSION_MAGIC)
    digest.update(f"scale={scale}\n".encode("ascii"))
    digest.update(f"rows={rows}\n".encode("ascii"))
    digest.update(f"columns={columns}\n".encode("ascii"))
    digest.update(b"column_ids\n")
    return digest


def _redacted_report(
    contract: dict[str, Any],
    observations: dict[str, Any],
    failures: Iterable[str],
) -> dict[str, Any]:
    failure_list = sorted(set(failures))
    required_scale = contract["semantic_digests"]["required_scale"]
    reference = contract["semantic_digests"]["reference_sha256"]
    semantic_hashes = observations.get("semantic_sha256", {})
    diagnostic_matches = {
        key: semantic_hashes.get(key) == expected
        for key, expected in reference.items()
        if key != f"fixed{required_scale}"
    }
    return {
        "semantic_contract_version": contract["semantic_contract_version"],
        "status": "pass" if not failure_list else "fail",
        "failed_checks": failure_list,
        "shape": {
            "rows": observations.get("rows"),
            "columns": observations.get("columns"),
            "values": observations.get("values"),
        },
        "ordered_identifier_sha256": {
            "rows": observations.get("row_ids_sha256"),
            "columns": observations.get("column_ids_sha256"),
        },
        "semantic_sha256": semantic_hashes,
        "required_compatibility": {
            "scale": required_scale,
            "match": semantic_hashes.get(f"fixed{required_scale}")
            == reference[f"fixed{required_scale}"],
        },
        "diagnostic_reference_matches": diagnostic_matches,
    }


def _write_json_atomic(path: Path, report: dict[str, Any]) -> None:
    if path.exists() and path.is_dir():
        raise VerificationError("Report path must name a file")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            json.dump(report, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def scan_expression(path: Path, contract: dict[str, Any]) -> dict[str, Any]:
    """Return aggregate observations while retaining no expression matrix."""

    if not path.is_file():
        raise VerificationError("Expression input is missing")
    expected_rows = contract["shape"]["rows"]
    expected_columns = contract["shape"]["columns"]
    required_scale = contract["semantic_digests"]["required_scale"]
    scales = [required_scale, *contract["semantic_digests"]["diagnostic_scales"]]
    digests = {
        scale: _initial_semantic_digest(scale, expected_rows, expected_columns)
        for scale in scales
    }
    row_digest = _sha256()
    column_digest = _sha256()
    row_digest.update(ORDERED_IDS_MAGIC)
    column_digest.update(ORDERED_IDS_MAGIC)
    row_ids: set[str] = set()
    rows = 0
    values = 0
    zero_count = 0
    minimum: Decimal | None = None
    maximum: Decimal | None = None

    try:
        handle = path.open("r", encoding="utf-8-sig", newline="")
    except (OSError, UnicodeError) as exc:
        raise VerificationError("Could not open expression input") from exc

    try:
        reader = csv.reader(handle, strict=True)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise VerificationError("Expression CSV is empty") from exc
        except csv.Error as exc:
            raise VerificationError("Expression CSV header is malformed") from exc
        if len(header) != expected_columns + 1:
            raise VerificationError("Expression CSV column count is incompatible")
        column_ids = header[1:]
        if any(not identifier for identifier in column_ids):
            raise VerificationError("Expression CSV contains an empty column identifier")
        if len(set(column_ids)) != len(column_ids):
            raise VerificationError("Expression CSV column identifiers are not unique")
        for identifier in column_ids:
            record = _identifier_record(identifier)
            column_digest.update(record)
            for digest in digests.values():
                digest.update(record)
        for digest in digests.values():
            digest.update(b"rows\n")

        with localcontext() as context:
            context.prec = 50
            quantums = {scale: Decimal(1).scaleb(-scale) for scale in scales}
            multipliers = {scale: Decimal(10) ** scale for scale in scales}
            try:
                for row in reader:
                    rows += 1
                    if rows > expected_rows:
                        raise VerificationError("Expression CSV has too many data rows")
                    if len(row) != expected_columns + 1:
                        raise VerificationError(
                            "Expression CSV has an incompatible data-row width"
                        )
                    row_id = row[0]
                    if not row_id:
                        raise VerificationError("Expression CSV contains an empty row identifier")
                    if row_id in row_ids:
                        raise VerificationError("Expression CSV row identifiers are not unique")
                    row_ids.add(row_id)
                    row_digest.update(_identifier_record(row_id))
                    prefixes = {scale: _row_prefix(row_id) for scale in scales}
                    scaled_fields: dict[int, list[bytes]] = {
                        scale: [] for scale in scales
                    }

                    for token in row[1:]:
                        values += 1
                        try:
                            number = Decimal(token.strip())
                        except (DecimalException, ValueError) as exc:
                            raise VerificationError(
                                "Expression CSV contains an invalid numeric field"
                            ) from exc
                        if not number.is_finite():
                            raise VerificationError(
                                "Expression CSV contains a non-finite numeric field"
                            )
                        if number < 0:
                            raise VerificationError(
                                "Expression CSV contains a negative numeric field"
                            )
                        if number.is_zero():
                            zero_count += 1
                        minimum = number if minimum is None or number < minimum else minimum
                        maximum = number if maximum is None or number > maximum else maximum
                        for scale in scales:
                            try:
                                scaled = int(
                                    number.quantize(
                                        quantums[scale], rounding=ROUND_HALF_UP
                                    )
                                    * multipliers[scale]
                                )
                            except (
                                DecimalException,
                                OverflowError,
                                ValueError,
                            ) as exc:
                                raise VerificationError(
                                    "Expression CSV contains a numeric field "
                                    "outside the supported decimal range"
                                ) from exc
                            scaled_fields[scale].append(
                                str(0 if scaled == 0 else scaled).encode("ascii")
                            )

                    for scale, digest in digests.items():
                        digest.update(prefixes[scale])
                        digest.update(b"\t")
                        digest.update(b"\t".join(scaled_fields[scale]))
                        digest.update(b"\n")
            except csv.Error as exc:
                raise VerificationError("Expression CSV data are malformed") from exc
    except UnicodeError as exc:
        raise VerificationError("Expression CSV is not valid UTF-8") from exc
    finally:
        handle.close()

    return {
        "rows": rows,
        "columns": len(column_ids),
        "values": values,
        "zero_count": zero_count,
        "minimum": str(minimum) if minimum is not None else None,
        "maximum": str(maximum) if maximum is not None else None,
        "row_ids_sha256": row_digest.hexdigest(),
        "column_ids_sha256": column_digest.hexdigest(),
        "semantic_sha256": {
            f"fixed{scale}": digest.hexdigest()
            for scale, digest in sorted(digests.items())
        },
    }


def verify_expression(
    path: Path,
    contract_path: Path = DEFAULT_CONTRACT,
    report_path: Path | None = None,
) -> dict[str, Any]:
    contract = load_contract(contract_path)
    try:
        observations = scan_expression(path, contract)
    except VerificationError:
        if report_path is not None:
            _write_json_atomic(
                report_path,
                {
                    "semantic_contract_version": contract["semantic_contract_version"],
                    "status": "fail",
                    "failed_checks": ["input_structure_or_numeric_validity"],
                    "shape": {"rows": None, "columns": None, "values": None},
                    "ordered_identifier_sha256": {
                        "rows": None,
                        "columns": None,
                    },
                    "semantic_sha256": {},
                    "required_compatibility": {
                        "scale": contract["semantic_digests"]["required_scale"],
                        "match": False,
                    },
                    "diagnostic_reference_matches": {},
                },
            )
        raise
    failures: list[str] = []
    shape = contract["shape"]
    identifiers = contract["ordered_identifier_digests"]
    numeric = contract["numeric_constraints"]
    semantics = contract["semantic_digests"]
    if observations["rows"] != shape["rows"]:
        failures.append("row_count")
    if observations["columns"] != shape["columns"]:
        failures.append("column_count")
    if observations["values"] != shape["values"]:
        failures.append("value_count")
    if observations["zero_count"] != numeric["zero_count"]:
        failures.append("zero_count")
    if observations["row_ids_sha256"] != identifiers["row_ids_sha256"]:
        failures.append("ordered_row_identifiers")
    if observations["column_ids_sha256"] != identifiers["column_ids_sha256"]:
        failures.append("ordered_column_identifiers")
    required_key = f"fixed{semantics['required_scale']}"
    if (
        observations["semantic_sha256"][required_key]
        != semantics["reference_sha256"][required_key]
    ):
        failures.append("required_fixed_scale_semantics")

    report = _redacted_report(contract, observations, failures)
    if report_path is not None:
        _write_json_atomic(report_path, report)
    if failures:
        raise VerificationError(
            "IMvigor210 expression semantic verification failed: "
            + ", ".join(sorted(failures))
        )
    return report


def run_self_test() -> None:
    # Non-sensitive fixture.  Constants were generated independently from the
    # byte framing documented in this module and guard accidental format drift.
    fixture = '"","sample-A","sample-B"\n"feature-1",0,1.2345675\n"feature-2",2.5,3\n'
    contract = {
        "semantic_contract_version": 1,
        "shape": {"rows": 2, "columns": 2, "values": 4},
        "numeric_constraints": {"zero_count": 1},
        "ordered_identifier_digests": {
            "row_ids_sha256": "55bddc56639c38463102a2639e29568893713aacbc976365e9b5b3ebeb713ee3",
            "column_ids_sha256": "d62695cafab034355274886969cef41f802981dfddbb1c3722c1051aeb6d263f",
        },
        "semantic_digests": {
            "required_scale": 6,
            "diagnostic_scales": [7, 8],
            "reference_sha256": {
                "fixed6": "3735d8755f800618217d0d930a6154ffec2c5091a1d75ae7ea4ca8b173c7367f",
                "fixed7": "79b7005704af91a0636bd0987c644f3ad2a19db79a9a1baeaaf19ea0bbd346dc",
                "fixed8": "be68a8454559a6e3175668db0932fe12f23190e0ba729e998f46bfbf7485cc3f",
            },
        },
    }
    with tempfile.TemporaryDirectory(prefix="imvigor210-semantic-self-test-") as temp:
        path = Path(temp) / "fixture.csv"
        with path.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write(fixture)
        observed = scan_expression(path, contract)
    expected = {
        "row_ids_sha256": "55bddc56639c38463102a2639e29568893713aacbc976365e9b5b3ebeb713ee3",
        "column_ids_sha256": "d62695cafab034355274886969cef41f802981dfddbb1c3722c1051aeb6d263f",
        "fixed6": "3735d8755f800618217d0d930a6154ffec2c5091a1d75ae7ea4ca8b173c7367f",
        "fixed7": "79b7005704af91a0636bd0987c644f3ad2a19db79a9a1baeaaf19ea0bbd346dc",
        "fixed8": "be68a8454559a6e3175668db0932fe12f23190e0ba729e998f46bfbf7485cc3f",
    }
    actual = {
        "row_ids_sha256": observed["row_ids_sha256"],
        "column_ids_sha256": observed["column_ids_sha256"],
        **observed["semantic_sha256"],
    }
    if actual != expected:
        raise VerificationError("Embedded semantic-framing self-test failed")


def main() -> int:
    arguments = parse_args()
    try:
        if arguments.self_test:
            run_self_test()
            print("IMvigor210 semantic framing self-test: pass")
            return 0
        report = verify_expression(
            arguments.input, arguments.contract, arguments.report
        )
        print(json.dumps(report, sort_keys=True, separators=(",", ":")))
        return 0
    except VerificationError as exc:
        print(f"ERROR: {exc}", file=os.sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
