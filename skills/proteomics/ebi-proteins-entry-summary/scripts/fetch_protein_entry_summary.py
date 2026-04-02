#!/usr/bin/env python3
"""Fetch a compact EBI Proteins API accession summary."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen


API_ROOT = "https://www.ebi.ac.uk/proteins/api/proteins"
DEFAULT_ACCESSION = "P38398"
USER_AGENT = "SciSkillUniverse/0.2"


def normalize_accession(accession: str) -> str:
    normalized = accession.strip().upper()
    if not normalized:
        raise ValueError("Accession must not be empty.")
    return normalized


def fetch_entry(accession: str) -> dict[str, Any]:
    normalized = normalize_accession(accession)
    request = Request(
        f"{API_ROOT}/{quote(normalized, safe='')}",
        headers={
            "Accept": "application/json",
            "User-Agent": USER_AGENT,
        },
    )
    try:
        with urlopen(request, timeout=30) as response:
            return json.load(response)
    except HTTPError as exc:
        message = exc.read().decode("utf-8", errors="replace").strip()
        raise SystemExit(f"EBI Proteins request failed for {normalized}: HTTP {exc.code} {message}") from exc
    except URLError as exc:
        raise SystemExit(f"EBI Proteins request failed for {normalized}: {exc.reason}") from exc


def _name_values(records: Any) -> list[str]:
    if not isinstance(records, list):
        return []
    values: list[str] = []
    for record in records:
        if not isinstance(record, dict):
            continue
        value = record.get("value")
        if isinstance(value, str) and value.strip() and value not in values:
            values.append(value.strip())
    return values


def _gene_names(records: Any) -> list[str]:
    if not isinstance(records, list):
        return []
    values: list[str] = []
    for record in records:
        if not isinstance(record, dict):
            continue
        primary = record.get("name")
        if isinstance(primary, dict):
            for value in _name_values([primary]):
                if value not in values:
                    values.append(value)
        for value in _name_values(record.get("synonyms")):
            if value not in values:
                values.append(value)
    return values


def _recommended_name(record: Any) -> str | None:
    if not isinstance(record, dict):
        return None
    recommended = record.get("recommendedName")
    if not isinstance(recommended, dict):
        return None
    full_name = recommended.get("fullName")
    if not isinstance(full_name, dict):
        return None
    value = full_name.get("value")
    if isinstance(value, str) and value.strip():
        return value.strip()
    return None


def _comment_entries(records: Any, limit: int = 3) -> list[dict[str, Any]]:
    if not isinstance(records, list):
        return []
    comments: list[dict[str, Any]] = []
    for record in records:
        if not isinstance(record, dict):
            continue
        text_blocks = record.get("text")
        text_value = None
        if isinstance(text_blocks, list):
            for block in text_blocks:
                if not isinstance(block, dict):
                    continue
                value = block.get("value")
                if isinstance(value, str) and value.strip():
                    text_value = value.strip()
                    break
        comments.append(
            {
                "type": record.get("type"),
                "text": text_value,
            }
        )
        if len(comments) >= limit:
            break
    return comments


def _feature_entries(records: Any, limit: int = 5) -> list[dict[str, Any]]:
    if not isinstance(records, list):
        return []
    features: list[dict[str, Any]] = []
    for record in records:
        if not isinstance(record, dict):
            continue
        features.append(
            {
                "type": record.get("type"),
                "category": record.get("category"),
                "description": record.get("description"),
                "begin": record.get("begin"),
                "end": record.get("end"),
            }
        )
        if len(features) >= limit:
            break
    return features


def build_summary(payload: dict[str, Any]) -> dict[str, Any]:
    organism = payload.get("organism") if isinstance(payload.get("organism"), dict) else {}
    names = organism.get("names") if isinstance(organism.get("names"), list) else []
    scientific_name = None
    common_name = None
    for record in names:
        if not isinstance(record, dict):
            continue
        name_type = record.get("type")
        value = record.get("value")
        if not isinstance(value, str) or not value.strip():
            continue
        if name_type == "scientific" and scientific_name is None:
            scientific_name = value.strip()
        if name_type == "common" and common_name is None:
            common_name = value.strip()

    sequence = payload.get("sequence") if isinstance(payload.get("sequence"), dict) else {}
    keywords = _name_values(payload.get("keywords"))[:10]
    accession = payload.get("accession") or DEFAULT_ACCESSION
    return {
        "accession": accession,
        "entry_id": payload.get("id"),
        "recommended_name": _recommended_name(payload.get("protein")),
        "gene_names": _gene_names(payload.get("gene")),
        "organism": scientific_name,
        "organism_scientific_name": scientific_name,
        "organism_common_name": common_name,
        "taxonomy_id": organism.get("taxonomy"),
        "sequence_length": sequence.get("length"),
        "sequence_mass": sequence.get("mass"),
        "keywords": keywords,
        "comments": _comment_entries(payload.get("comments")),
        "features": _feature_entries(payload.get("features")),
        "source_url": f"{API_ROOT}/{quote(str(accession), safe='')}",
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
    parser.add_argument("--accession", default=DEFAULT_ACCESSION, help="Protein accession such as P38398.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    payload = fetch_entry(args.accession)
    summary = build_summary(payload)
    write_json(summary, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
