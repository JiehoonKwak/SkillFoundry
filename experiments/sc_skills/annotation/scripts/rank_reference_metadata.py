#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from reference_search import build_reference_selection_summary, rank_reference_candidates


def _read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Rank exported CELLxGENE/Census metadata for spatial annotation reference selection.")
    parser.add_argument("--metadata", type=Path, required=True, help="Input CSV/TSV with dataset metadata rows.")
    parser.add_argument("--species", required=True)
    parser.add_argument("--tissue", required=True)
    parser.add_argument("--condition", default="normal")
    parser.add_argument("--min-cells", type=int, default=1000)
    parser.add_argument("--out", type=Path, required=True, help="Ranked TSV output path.")
    parser.add_argument("--summary-out", type=Path, default=None, help="Optional JSON selection summary path.")
    args = parser.parse_args(argv)

    frame = _read_table(args.metadata)
    ranked = rank_reference_candidates(
        frame.to_dict(orient="records"),
        {
            "species": args.species,
            "tissue": args.tissue,
            "condition": args.condition,
            "min_cells": args.min_cells,
        },
    )
    summary = build_reference_selection_summary(
        ranked,
        {
            "species": args.species,
            "tissue": args.tissue,
            "condition": args.condition,
            "min_cells": args.min_cells,
        },
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    ranked.to_csv(args.out, sep="\t", index=False)
    if args.summary_out is not None:
        args.summary_out.parent.mkdir(parents=True, exist_ok=True)
        args.summary_out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(
        json.dumps(
            {
                "candidate_count": int(ranked.shape[0]),
                "accepted_candidate_count": int(ranked["accept_for_transfer"].sum()),
                "selected_dataset_id": summary["selected_dataset_id"],
                "out": str(args.out.resolve()),
                "summary_out": str(args.summary_out.resolve()) if args.summary_out is not None else None,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
