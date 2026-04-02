#!/usr/bin/env python3
"""Fetch an RCSB PDB entry summary by accession."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen


API_ROOT = "https://data.rcsb.org/rest/v1/core/entry"


def fetch_json(entry_id: str) -> dict:
    request = Request(f"{API_ROOT}/{quote(entry_id.upper())}", headers={"User-Agent": "SciSkillUniverse/0.2"})
    with urlopen(request, timeout=30) as response:
        return json.load(response)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--entry-id", required=True, help="PDB accession, for example 4HHB")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path")
    args = parser.parse_args()

    payload = fetch_json(args.entry_id)
    text = json.dumps(payload, indent=2, sort_keys=True)
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
