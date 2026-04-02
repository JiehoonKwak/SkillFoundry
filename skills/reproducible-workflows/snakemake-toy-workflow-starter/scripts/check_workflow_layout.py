#!/usr/bin/env python3
"""Validate the toy Snakemake starter layout."""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"


def main() -> int:
    required = [
        EXAMPLES / "Snakefile",
        EXAMPLES / "config.yaml",
        EXAMPLES / "data" / "input.txt",
    ]
    missing = [path for path in required if not path.exists()]
    if missing:
        for path in missing:
            print(f"missing: {path}")
        return 1
    print("Snakemake toy workflow starter layout looks complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
