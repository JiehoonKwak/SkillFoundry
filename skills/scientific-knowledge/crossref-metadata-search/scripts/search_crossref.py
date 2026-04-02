#!/usr/bin/env python3
"""Search Crossref works metadata by title query."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen


API_URL = "https://api.crossref.org/works"


def build_url(query_title: str, rows: int) -> str:
    return f"{API_URL}?{urlencode({'query.title': query_title, 'rows': rows})}"


def fetch_json(url: str) -> dict:
    request = Request(url, headers={"User-Agent": "SciSkillUniverse/0.1 (mailto:test@example.com)"})
    with urlopen(request, timeout=30) as response:
        return json.load(response)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query-title", required=True, help="Crossref title query")
    parser.add_argument("--rows", type=int, default=5, help="Number of rows to request")
    parser.add_argument("--out", type=Path, default=None, help="Optional output path")
    args = parser.parse_args()

    if args.rows < 1 or args.rows > 100:
        raise SystemExit("--rows must be between 1 and 100")

    payload = fetch_json(build_url(args.query_title, args.rows))
    text = json.dumps(payload, indent=2, sort_keys=True)
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
