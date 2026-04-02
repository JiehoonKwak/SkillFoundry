#!/usr/bin/env python3
"""Search the current ClinicalTrials.gov v2 API."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen


API_URL = "https://clinicaltrials.gov/api/v2/studies"


def build_url(condition: str | None, term: str | None, page_size: int) -> str:
    params = {"pageSize": page_size}
    if condition:
        params["query.cond"] = condition
    if term:
        params["query.term"] = term
    return f"{API_URL}?{urlencode(params)}"


def fetch_json(url: str) -> dict:
    request = Request(url, headers={"User-Agent": "SciSkillUniverse/0.2"})
    with urlopen(request, timeout=30) as response:
        return json.load(response)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--condition", help="Condition query, for example melanoma")
    group.add_argument("--term", help="General term query")
    parser.add_argument("--page-size", type=int, default=5, help="Maximum number of studies to request")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path")
    args = parser.parse_args()

    if args.page_size < 1 or args.page_size > 20:
        raise SystemExit("--page-size must be between 1 and 20")

    payload = fetch_json(build_url(args.condition, args.term, args.page_size))
    text = json.dumps(payload, indent=2, sort_keys=True)
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
