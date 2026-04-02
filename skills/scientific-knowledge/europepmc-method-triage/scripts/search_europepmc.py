#!/usr/bin/env python3
"""Small Europe PMC search helper for method triage."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen


API_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


def build_url(query: str, page_size: int, sort: str | None) -> str:
    params = {"query": query, "pageSize": page_size, "format": "json"}
    if sort:
        params["sort"] = sort
    return f"{API_URL}?{urlencode(params)}"


def fetch_json(url: str) -> dict:
    request = Request(url, headers={"User-Agent": "SciSkillUniverse/0.1"})
    with urlopen(request, timeout=30) as response:
        return json.load(response)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query", required=True, help="Europe PMC query string")
    parser.add_argument("--page-size", type=int, default=5, help="Number of results to request")
    parser.add_argument("--sort", default=None, help="Optional Europe PMC sort value")
    parser.add_argument("--out", type=Path, default=None, help="Optional path to save the JSON payload")
    args = parser.parse_args()

    if args.page_size < 1 or args.page_size > 1000:
        raise SystemExit("--page-size must be between 1 and 1000")

    payload = fetch_json(build_url(args.query, args.page_size, args.sort))
    text = json.dumps(payload, indent=2, sort_keys=True)
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
