#!/usr/bin/env python3
"""Download checksum-pinned public inputs into the ignored local cache."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import ssl
import sys
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "data" / "source_manifest.tsv"
DEFAULT_DESTINATION = ROOT / "data" / "raw"
BUFFER_SIZE = 1024 * 1024
USER_AGENT = "car-t-bystander-resistance-reproducibility/1.0"


@dataclass(frozen=True)
class Source:
    source_id: str
    cohort_id: str
    citation: str
    source_url: str
    expected_filename: str
    sha256: str
    size_bytes: int
    access_class: str
    license_name: str
    license_url: str
    fetch_mode: str
    repository_policy: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fetch public processed sources listed in data/source_manifest.tsv. "
            "Files are written only to the ignored local cache."
        )
    )
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--destination", type=Path, default=DEFAULT_DESTINATION)
    parser.add_argument(
        "--source",
        action="append",
        default=[],
        help="Source ID to fetch; repeat for multiple sources. Default: automatic rows.",
    )
    parser.add_argument("--list", action="store_true", help="List manifest sources and exit.")
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="Verify existing files without downloading missing inputs.",
    )
    parser.add_argument(
        "--accept-licensed-public-downloads",
        action="store_true",
        help=(
            "Permit explicitly selected manual_public rows after reviewing their "
            "licence and attribution requirements."
        ),
    )
    parser.add_argument("--force", action="store_true", help="Replace a verified file.")
    return parser.parse_args()


def read_manifest(path: Path) -> list[Source]:
    if not path.is_file():
        raise FileNotFoundError(f"Source manifest not found: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "source_id",
            "cohort_id",
            "citation",
            "source_url",
            "expected_filename",
            "sha256",
            "size_bytes",
            "access_class",
            "license",
            "license_url",
            "fetch_mode",
            "repository_policy",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Manifest lacks columns: {', '.join(sorted(missing))}")
        rows = list(reader)

    sources: list[Source] = []
    seen: set[str] = set()
    for line_number, row in enumerate(rows, start=2):
        source_id = row["source_id"].strip()
        if not source_id or source_id in seen:
            raise ValueError(f"Invalid or duplicate source_id on manifest line {line_number}")
        seen.add(source_id)
        filename = row["expected_filename"].strip()
        if not filename or Path(filename).name != filename:
            raise ValueError(f"Unsafe expected_filename for {source_id}: {filename!r}")
        digest = row["sha256"].strip().lower()
        if len(digest) != 64 or any(ch not in "0123456789abcdef" for ch in digest):
            raise ValueError(f"Invalid SHA-256 for {source_id}")
        try:
            size_bytes = int(row["size_bytes"])
        except ValueError as exc:
            raise ValueError(f"Invalid size_bytes for {source_id}") from exc
        sources.append(
            Source(
                source_id=source_id,
                cohort_id=row["cohort_id"].strip(),
                citation=row["citation"].strip(),
                source_url=row["source_url"].strip(),
                expected_filename=filename,
                sha256=digest,
                size_bytes=size_bytes,
                access_class=row["access_class"].strip(),
                license_name=row["license"].strip(),
                license_url=row["license_url"].strip(),
                fetch_mode=row["fetch_mode"].strip(),
                repository_policy=row["repository_policy"].strip(),
            )
        )
    return sources


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(BUFFER_SIZE), b""):
            digest.update(block)
    return digest.hexdigest()


def verify(path: Path, source: Source) -> tuple[bool, str]:
    if not path.is_file():
        return False, "missing"
    observed_size = path.stat().st_size
    if observed_size != source.size_bytes:
        return False, f"size mismatch ({observed_size} != {source.size_bytes})"
    observed_digest = file_sha256(path)
    if observed_digest != source.sha256:
        return False, f"SHA-256 mismatch ({observed_digest})"
    return True, "verified"


def download(source: Source, destination: Path, force: bool) -> Path:
    if not source.source_url:
        raise ValueError(f"{source.source_id} is a local input and has no download URL")
    destination.mkdir(parents=True, exist_ok=True)
    final_path = destination / source.expected_filename
    ok, status = verify(final_path, source)
    if ok and not force:
        print(f"VERIFIED\t{source.source_id}\t{final_path}")
        return final_path
    if final_path.exists() and not force:
        raise RuntimeError(
            f"Refusing to overwrite unverified file {final_path}: {status}. "
            "Remove it or use --force after inspection."
        )

    part_path = final_path.with_name(final_path.name + ".part")
    if part_path.exists():
        part_path.unlink()
    request = urllib.request.Request(source.source_url, headers={"User-Agent": USER_AGENT})
    context = ssl.create_default_context()
    digest = hashlib.sha256()
    bytes_written = 0
    try:
        with urllib.request.urlopen(request, timeout=120, context=context) as response:
            if response.status != 200:
                raise RuntimeError(f"HTTP {response.status} for {source.source_url}")
            with part_path.open("xb") as handle:
                while True:
                    block = response.read(BUFFER_SIZE)
                    if not block:
                        break
                    handle.write(block)
                    digest.update(block)
                    bytes_written += len(block)
                handle.flush()
                os.fsync(handle.fileno())
        if bytes_written != source.size_bytes:
            raise RuntimeError(
                f"Size mismatch for {source.source_id}: {bytes_written} != {source.size_bytes}"
            )
        observed_digest = digest.hexdigest()
        if observed_digest != source.sha256:
            raise RuntimeError(
                f"SHA-256 mismatch for {source.source_id}: {observed_digest}"
            )
        os.replace(part_path, final_path)
    except BaseException:
        if part_path.exists():
            part_path.unlink()
        raise
    print(f"DOWNLOADED\t{source.source_id}\t{final_path}")
    return final_path


def print_sources(sources: list[Source]) -> None:
    print("source_id\tcohort_id\tfetch_mode\tfilename\tlicence\trepository_policy")
    for source in sources:
        print(
            "\t".join(
                [
                    source.source_id,
                    source.cohort_id,
                    source.fetch_mode,
                    source.expected_filename,
                    source.license_name,
                    source.repository_policy,
                ]
            )
        )


def select_sources(sources: list[Source], args: argparse.Namespace) -> list[Source]:
    by_id = {source.source_id: source for source in sources}
    if args.source:
        unknown = sorted(set(args.source).difference(by_id))
        if unknown:
            raise ValueError(f"Unknown source ID(s): {', '.join(unknown)}")
        selected = [by_id[source_id] for source_id in dict.fromkeys(args.source)]
    else:
        selected = [source for source in sources if source.fetch_mode == "automatic"]

    for source in selected:
        if source.fetch_mode == "local_input":
            raise ValueError(
                f"{source.source_id} is a project-derived local input; follow data/README.md"
            )
        if source.fetch_mode == "manual_public" and not args.accept_licensed_public_downloads:
            raise ValueError(
                f"{source.source_id} requires --accept-licensed-public-downloads after "
                f"reviewing {source.license_name}: {source.license_url}"
            )
    return selected


def main() -> int:
    args = parse_args()
    try:
        sources = read_manifest(args.manifest.resolve())
        if args.list:
            print_sources(sources)
            return 0
        selected = select_sources(sources, args)
        if not selected:
            raise ValueError("No downloadable sources selected")
        destination = args.destination.resolve()
        for source in selected:
            target = destination / source.expected_filename
            if args.verify_only:
                ok, status = verify(target, source)
                label = "VERIFIED" if ok else "FAILED"
                print(f"{label}\t{source.source_id}\t{target}\t{status}")
                if not ok:
                    return 1
            else:
                download(source, destination, args.force)
        return 0
    except (OSError, ValueError, RuntimeError, urllib.error.URLError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
