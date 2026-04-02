#!/usr/bin/env python3
"""Search RCSB PDB entries with a minimal full-text query."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from urllib.request import Request, urlopen


API_URL = "https://search.rcsb.org/rcsbsearch/v2/query"


def build_payload(query: str, rows: int) -> bytes:
    payload = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {"value": query},
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": rows}},
    }
    return json.dumps(payload).encode("utf-8")


def fetch_json(payload: bytes) -> dict:
    request = Request(
        API_URL,
        data=payload,
        headers={"Content-Type": "application/json", "User-Agent": "SciSkillUniverse/0.1"},
    )
    with urlopen(request, timeout=30) as response:
        return json.load(response)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query", required=True, help="RCSB full-text query")
    parser.add_argument("--rows", type=int, default=5, help="Number of rows to request")
    parser.add_argument("--out", type=Path, default=None, help="Optional output path")
    args = parser.parse_args()

    if args.rows < 1 or args.rows > 100:
        raise SystemExit("--rows must be between 1 and 100")

    payload = fetch_json(build_payload(args.query, args.rows))
    text = json.dumps(payload, indent=2, sort_keys=True)
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
