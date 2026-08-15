#!/usr/bin/env python3
"""Record runtime and numerical-library provenance for an RNA-seq CI run."""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import os
import platform
import sys
import zlib
from pathlib import Path

import numpy as np


THREAD_ENVIRONMENT_VARIABLES = (
    "OMP_NUM_THREADS",
    "OMP_DYNAMIC",
    "OMP_PROC_BIND",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "BLIS_NUM_THREADS",
)


def numpy_configuration() -> str:
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        np.show_config()
    return buffer.getvalue().rstrip()


def collect_provenance() -> dict[str, object]:
    return {
        "schema_version": 1,
        "platform": {
            "platform": platform.platform(),
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
            "processor": platform.processor(),
            "cpu_count": os.cpu_count(),
        },
        "python": {
            "implementation": platform.python_implementation(),
            "version": platform.python_version(),
            "version_string": sys.version,
        },
        "zlib": {
            "compiled_version": zlib.ZLIB_VERSION,
            "runtime_version": zlib.ZLIB_RUNTIME_VERSION,
        },
        "numpy": {
            "version": np.__version__,
            "show_config": numpy_configuration(),
        },
        "thread_environment": {
            name: os.environ.get(name) for name in THREAD_ENVIRONMENT_VARIABLES
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(collect_provenance(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
