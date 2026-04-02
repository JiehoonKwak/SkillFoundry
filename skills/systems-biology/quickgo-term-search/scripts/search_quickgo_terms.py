#!/usr/bin/env python3
"""Search Gene Ontology terms through the official QuickGO API."""

from __future__ import annotations

import argparse
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any


API_ROOT = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/search"
USER_AGENT = "SciSkillUniverse/1.0 (+https://example.invalid)"


def normalize_positive_int(value: int, label: str) -> int:
    if value <= 0:
        raise ValueError(f"{label} must be a positive integer.")
    return value


def compact_term_record(record: dict[str, Any]) -> dict[str, Any]:
    definition = record.get("definition") or {}
    return {
        "id": record.get("id"),
        "name": record.get("name"),
        "aspect": record.get("aspect"),
        "is_obsolete": bool(record.get("isObsolete", False)),
        "definition": definition.get("text"),
    }


def fetch_json(url: str, attempts: int = 3, timeout: int = 30) -> dict[str, Any]:
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        request = urllib.request.Request(
            url,
            headers={
                "Accept": "application/json",
                "User-Agent": USER_AGENT,
            },
        )
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return json.load(response)
        except urllib.error.HTTPError as exc:
            last_error = exc
            if exc.code not in {429, 500, 502, 503, 504} or attempt == attempts:
                raise
        except urllib.error.URLError as exc:
            last_error = exc
            if attempt == attempts:
                raise
        time.sleep(min(2 ** attempt, 5))
    if last_error is not None:
        raise last_error
    raise RuntimeError("QuickGO request failed without a captured error.")


def search_quickgo_terms(query: str, limit: int) -> dict[str, Any]:
    params = urllib.parse.urlencode({"query": query, "page": 1, "limit": limit})
    payload = fetch_json(f"{API_ROOT}?{params}")
    results = payload.get("results") or []
    terms = [compact_term_record(record) for record in results[:limit]]
    return {
        "query": query,
        "term_count": len(terms),
        "number_of_hits": payload.get("numberOfHits", 0),
        "terms": terms,
    }


def write_json(payload: dict[str, Any], out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query", required=True, help="Free-text GO term query.")
    parser.add_argument("--limit", type=int, default=5, help="Maximum number of terms to return.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    try:
        limit = normalize_positive_int(args.limit, "limit")
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    payload = search_quickgo_terms(query=args.query, limit=limit)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
