#!/usr/bin/env python3
"""Search the CELLxGENE Census dataset catalog by keyword."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load_stack():
    try:
        import cellxgene_census
        import pandas as pd
    except ImportError as exc:  # pragma: no cover - exercised through runtime tests
        raise SystemExit("This script requires the census prefix at slurm/envs/census.") from exc
    return cellxgene_census, pd


def search_census(keyword: str, limit: int, census_version: str) -> dict:
    cellxgene_census, _ = load_stack()
    keyword = keyword.strip()
    if not keyword:
        raise SystemExit("keyword must not be empty.")
    census = cellxgene_census.open_soma(census_version=census_version)
    try:
        datasets = census["census_info"]["datasets"].read().concat().to_pandas()
    finally:
        census.close()

    mask = datasets["dataset_title"].str.contains(keyword, case=False, na=False) | datasets["collection_name"].str.contains(
        keyword, case=False, na=False
    )
    matches = datasets.loc[
        mask,
        ["dataset_id", "dataset_title", "collection_name", "dataset_total_cell_count"],
    ].copy()
    if matches.empty:
        raise SystemExit(f"No CELLxGENE Census datasets matched keyword: {keyword}")
    matches = matches.sort_values("dataset_total_cell_count", ascending=False).head(limit)
    top_matches = [
        {
            "dataset_id": row["dataset_id"],
            "dataset_title": row["dataset_title"],
            "collection_name": row["collection_name"],
            "dataset_total_cell_count": int(row["dataset_total_cell_count"]),
        }
        for _, row in matches.iterrows()
    ]
    return {
        "keyword": keyword,
        "census_version": census_version,
        "match_count": int(mask.sum()),
        "top_match_count": len(top_matches),
        "unique_collections_in_top_matches": int(matches["collection_name"].nunique()),
        "total_cells_across_top_matches": int(matches["dataset_total_cell_count"].sum()),
        "top_matches": top_matches,
    }


def write_json(payload: dict, out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--keyword", default="Tabula Sapiens", help="Keyword for dataset_title or collection_name.")
    parser.add_argument("--limit", type=int, default=5, help="Number of top matches to keep.")
    parser.add_argument("--census-version", default="2025-11-08", help="Pinned Census version.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()
    if args.limit <= 0:
        raise SystemExit("limit must be a positive integer.")
    payload = search_census(keyword=args.keyword, limit=args.limit, census_version=args.census_version)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
