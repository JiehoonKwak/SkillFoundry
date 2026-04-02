#!/usr/bin/env python3
"""Search GBIF dataset metadata and summarize a small result set."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


API_URL = "https://api.gbif.org/v1/dataset/search"
USER_AGENT = "SciSkillUniverse/0.2"
CANONICAL_QUERY = ("puma", 3)
ASSET_PATH = Path(__file__).resolve().parents[1] / "assets" / "puma_dataset_search.json"


def fetch_json(url: str, retries: int = 3) -> dict:
    for attempt in range(retries):
        request = Request(url, headers={"User-Agent": USER_AGENT})
        try:
            with urlopen(request, timeout=30) as response:
                return json.load(response)
        except (HTTPError, URLError) as exc:
            if attempt == retries - 1:
                raise exc
            time.sleep(1.0 + attempt)
    raise AssertionError("unreachable")  # pragma: no cover


def build_summary(query: str, limit: int) -> dict[str, object]:
    url = f"{API_URL}?{urlencode({'q': query, 'limit': limit})}"
    payload = fetch_json(url)
    results = payload.get("results", [])
    dataset_summaries = [
        {
            "dataset_key": record.get("key"),
            "title": record.get("title"),
            "doi": record.get("doi"),
            "license": record.get("license"),
            "publishing_country": record.get("publishingCountry"),
            "type": record.get("type"),
        }
        for record in results
    ]
    licenses_seen = sorted({item["license"] for item in dataset_summaries if item["license"]})
    publishing_countries = sorted(
        {item["publishing_country"] for item in dataset_summaries if item["publishing_country"]}
    )
    return {
        "query": query,
        "limit": int(limit),
        "result_count": len(dataset_summaries),
        "end_of_records": bool(payload.get("endOfRecords", False)),
        "first_dataset_key": dataset_summaries[0]["dataset_key"] if dataset_summaries else None,
        "first_title": dataset_summaries[0]["title"] if dataset_summaries else None,
        "licenses_seen": licenses_seen,
        "publishing_countries": publishing_countries,
        "dataset_summaries": dataset_summaries,
    }


def load_fallback() -> dict[str, object]:
    return json.loads(ASSET_PATH.read_text(encoding="utf-8"))


def build_summary_with_fallback(query: str, limit: int) -> dict[str, object]:
    try:
        return build_summary(query=query, limit=limit)
    except (HTTPError, URLError):
        if (query, limit) == CANONICAL_QUERY and ASSET_PATH.exists():
            payload = load_fallback()
            payload["fallback_used"] = True
            return payload
        raise


def write_json(payload: dict[str, object], out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query", required=True, help="GBIF dataset search query.")
    parser.add_argument("--limit", type=int, default=3, help="Maximum number of datasets to request.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    if args.limit < 1 or args.limit > 20:
        raise SystemExit("--limit must be between 1 and 20.")

    payload = build_summary_with_fallback(query=args.query, limit=args.limit)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
