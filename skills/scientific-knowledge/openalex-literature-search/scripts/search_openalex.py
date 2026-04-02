#!/usr/bin/env python3
"""Small OpenAlex search helper for literature triage."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen


API_URL = "https://api.openalex.org/works"


def build_url(query: str, per_page: int, mailto: str | None) -> str:
    params = {"search": query, "per-page": per_page}
    if mailto:
        params["mailto"] = mailto
    return f"{API_URL}?{urlencode(params)}"


def fetch_json(url: str) -> dict:
    request = Request(url, headers={"User-Agent": "SciSkillUniverse/0.1"})
    with urlopen(request, timeout=30) as response:
        return json.load(response)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query", required=True, help="OpenAlex search query")
    parser.add_argument("--per-page", type=int, default=5, help="Number of results to request")
    parser.add_argument("--mailto", default=None, help="Optional email for polite-pool requests")
    parser.add_argument("--out", type=Path, default=None, help="Optional path to save the JSON payload")
    args = parser.parse_args()

    if args.per_page < 1 or args.per_page > 200:
        raise SystemExit("--per-page must be between 1 and 200")

    payload = fetch_json(build_url(args.query, args.per_page, args.mailto))
    text = json.dumps(payload, indent=2, sort_keys=True)
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
