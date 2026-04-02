#!/usr/bin/env python3
"""Fetch a compact InterPro entry summary."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen


BASE_URL = "https://www.ebi.ac.uk/interpro/api/entry/interpro"


def fetch_entry(accession: str) -> dict:
    accession = accession.upper().strip()
    request = Request(
        f"{BASE_URL}/{accession}/",
        headers={"Accept": "application/json", "User-Agent": "SciSkillUniverse/1.0"},
    )
    try:
        with urlopen(request, timeout=60) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except HTTPError as exc:
        if exc.code == 404:
            raise SystemExit(f"InterPro accession not found: {accession}") from exc
        raise SystemExit(f"InterPro request failed with HTTP {exc.code}") from exc

    metadata = payload["metadata"]
    member_databases = metadata.get("member_databases") or {}
    go_terms = metadata.get("go_terms") or []
    name_block = metadata.get("name")
    if isinstance(name_block, dict):
        entry_name = name_block.get("name") or name_block.get("short") or accession
    else:
        entry_name = str(name_block or accession)
    hierarchy = metadata.get("hierarchy") or {}
    return {
        "accession": metadata["accession"],
        "name": entry_name,
        "type": metadata["type"],
        "source_database": metadata["source_database"],
        "member_database_count": len(member_databases),
        "member_database_names": sorted(member_databases.keys()),
        "go_term_count": len(go_terms),
        "go_term_identifiers": [term["identifier"] for term in go_terms[:5]],
        "hierarchy_name": hierarchy.get("name"),
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
    parser.add_argument("--accession", default="IPR000023", help="InterPro accession.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()
    payload = fetch_entry(args.accession)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
