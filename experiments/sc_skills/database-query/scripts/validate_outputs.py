#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from run_database_query import (
    DEFAULT_INPUT,
    load_json,
    read_tsv,
    validate_expected_invariants,
    validate_metadata_outputs,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Validate database-query starter outputs.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)

    metadata = load_json(Path(__file__).resolve().parents[1] / "metadata.yaml")
    validate_metadata_outputs(metadata, args.outdir)

    payload = load_json(args.input)
    results_rows = read_tsv(args.outdir / "query_results.tsv")
    normalization_rows = read_tsv(args.outdir / "normalization_qc.tsv")
    reconciliation_rows = read_tsv(args.outdir / "reconciliation_qc.tsv")
    run_summary = json.loads((args.outdir / "run_summary.json").read_text(encoding="utf-8"))
    validate_expected_invariants(payload, results_rows, normalization_rows, reconciliation_rows, run_summary)
    print(json.dumps({"validated": str(args.outdir.resolve())}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
