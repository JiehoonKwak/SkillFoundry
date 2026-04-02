#!/usr/bin/env python3
"""Fetch a compact UniProt sequence-feature summary."""

from __future__ import annotations

import argparse
import json
import time
from collections import Counter
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen


API_ROOT = "https://rest.uniprot.org/uniprotkb"
DEFAULT_ACCESSION = "P04637"
USER_AGENT = "SciSkillUniverse/1.0"
RETRYABLE_HTTP_CODES = {429, 500, 502, 503, 504}
FEATURE_PRIORITY = [
    "Chain",
    "Domain",
    "DNA binding",
    "Region",
    "Motif",
    "Binding site",
    "Site",
    "Modified residue",
    "Natural variant",
]


def normalize_accession(accession: str) -> str:
    normalized = accession.strip().upper()
    if not normalized:
        raise ValueError("Accession must not be empty.")
    return normalized


def fetch_entry(accession: str, attempts: int = 3) -> dict[str, Any]:
    normalized = normalize_accession(accession)
    request = Request(
        f"{API_ROOT}/{quote(normalized, safe='')}.json",
        headers={
            "Accept": "application/json",
            "User-Agent": USER_AGENT,
        },
    )
    for attempt in range(1, attempts + 1):
        try:
            with urlopen(request, timeout=60) as response:
                return json.load(response)
        except HTTPError as exc:
            message = exc.read().decode("utf-8", errors="replace").strip()
            if exc.code == 404:
                raise SystemExit(f"UniProt accession not found: {normalized}") from exc
            if exc.code == 400 and "invalid format" in message.lower():
                raise SystemExit(f"Invalid UniProt accession format: {normalized}") from exc
            if exc.code in RETRYABLE_HTTP_CODES and attempt < attempts:
                time.sleep(attempt)
                continue
            suffix = f" {message}" if message else ""
            raise SystemExit(f"UniProt request failed for {normalized}: HTTP {exc.code}{suffix}") from exc
        except URLError as exc:
            if attempt < attempts:
                time.sleep(attempt)
                continue
            raise SystemExit(f"UniProt request failed for {normalized}: {exc.reason}") from exc
    raise SystemExit(f"UniProt request failed for {normalized}: exhausted retries")


def _recommended_name(payload: dict[str, Any]) -> str | None:
    protein = payload.get("proteinDescription")
    if not isinstance(protein, dict):
        return None
    recommended = protein.get("recommendedName")
    if not isinstance(recommended, dict):
        return None
    full_name = recommended.get("fullName")
    if not isinstance(full_name, dict):
        return None
    value = full_name.get("value")
    if isinstance(value, str) and value.strip():
        return value.strip()
    return None


def _gene_names(payload: dict[str, Any]) -> list[str]:
    genes = payload.get("genes")
    if not isinstance(genes, list):
        return []
    values: list[str] = []
    for record in genes:
        if not isinstance(record, dict):
            continue
        gene_name = record.get("geneName")
        if isinstance(gene_name, dict):
            value = gene_name.get("value")
            if isinstance(value, str) and value.strip() and value.strip() not in values:
                values.append(value.strip())
        synonyms = record.get("synonyms")
        if isinstance(synonyms, list):
            for synonym in synonyms:
                if not isinstance(synonym, dict):
                    continue
                value = synonym.get("value")
                if isinstance(value, str) and value.strip() and value.strip() not in values:
                    values.append(value.strip())
    return values


def _location_point(payload: Any) -> tuple[int | None, str | None]:
    if not isinstance(payload, dict):
        return None, None
    value = payload.get("value")
    modifier = payload.get("modifier")
    if isinstance(value, int):
        return value, modifier if isinstance(modifier, str) else None
    return None, modifier if isinstance(modifier, str) else None


def _format_span(start: int | None, end: int | None) -> str | None:
    if start is None and end is None:
        return None
    if start is not None and end is not None:
        if start == end:
            return str(start)
        return f"{start}-{end}"
    if start is not None:
        return str(start)
    return str(end)


def _variation_label(feature: dict[str, Any]) -> str | None:
    alt = feature.get("alternativeSequence")
    if not isinstance(alt, dict):
        return None
    original = alt.get("originalSequence")
    alternatives = alt.get("alternativeSequences")
    if not isinstance(original, str) or not isinstance(alternatives, list):
        return None
    normalized_alternatives = [value for value in alternatives if isinstance(value, str) and value]
    if not normalized_alternatives:
        return None
    return f"{original}->{'/'.join(normalized_alternatives)}"


def normalize_feature(feature: dict[str, Any]) -> dict[str, Any]:
    location = feature.get("location") if isinstance(feature.get("location"), dict) else {}
    start, start_modifier = _location_point(location.get("start"))
    end, end_modifier = _location_point(location.get("end"))
    normalized: dict[str, Any] = {
        "type": feature.get("type") or "unknown",
        "description": (feature.get("description") or None) if isinstance(feature.get("description"), str) else None,
        "start": start,
        "end": end,
        "span": _format_span(start, end),
        "feature_id": feature.get("featureId") or None,
        "evidence_count": len(feature.get("evidences") or []),
    }
    if start_modifier and start_modifier != "EXACT":
        normalized["start_modifier"] = start_modifier
    if end_modifier and end_modifier != "EXACT":
        normalized["end_modifier"] = end_modifier
    variation = _variation_label(feature)
    if variation is not None:
        normalized["variation"] = variation
    return normalized


def select_representative_features(features: list[dict[str, Any]], limit: int = 8) -> list[dict[str, Any]]:
    normalized = [normalize_feature(feature) for feature in features if isinstance(feature, dict)]
    normalized.sort(
        key=lambda item: (
            item["start"] if item["start"] is not None else 10**9,
            item["end"] if item["end"] is not None else 10**9,
            item["type"],
            item.get("description") or "",
            item.get("feature_id") or "",
        )
    )
    type_counts = Counter(item["type"] for item in normalized)
    first_by_type: dict[str, dict[str, Any]] = {}
    for item in normalized:
        first_by_type.setdefault(item["type"], item)

    selected: list[dict[str, Any]] = []
    seen_types: set[str] = set()
    for feature_type in FEATURE_PRIORITY:
        if feature_type in first_by_type:
            selected.append(first_by_type[feature_type])
            seen_types.add(feature_type)
            if len(selected) >= limit:
                return selected

    remaining_types = sorted(first_by_type, key=lambda value: (-type_counts[value], value))
    for feature_type in remaining_types:
        if feature_type in seen_types:
            continue
        selected.append(first_by_type[feature_type])
        if len(selected) >= limit:
            break
    return selected


def build_summary(payload: dict[str, Any], feature_limit: int = 8) -> dict[str, Any]:
    accession = payload.get("primaryAccession") or DEFAULT_ACCESSION
    organism = payload.get("organism") if isinstance(payload.get("organism"), dict) else {}
    sequence = payload.get("sequence") if isinstance(payload.get("sequence"), dict) else {}
    features = payload.get("features") if isinstance(payload.get("features"), list) else []
    feature_counter = Counter(
        feature.get("type") or "unknown"
        for feature in features
        if isinstance(feature, dict)
    )
    top_feature_types = [
        {"type": feature_type, "count": count}
        for feature_type, count in sorted(feature_counter.items(), key=lambda item: (-item[1], item[0]))[:10]
    ]
    return {
        "accession": accession,
        "entry_id": payload.get("uniProtkbId"),
        "recommended_name": _recommended_name(payload),
        "gene_names": _gene_names(payload),
        "organism_scientific_name": organism.get("scientificName"),
        "organism_common_name": organism.get("commonName"),
        "taxonomy_id": organism.get("taxonId"),
        "sequence_length": sequence.get("length"),
        "annotation_score": payload.get("annotationScore"),
        "feature_count": len(features),
        "feature_type_counts": dict(sorted(feature_counter.items())),
        "top_feature_types": top_feature_types,
        "representative_features": select_representative_features(features, limit=feature_limit),
        "source_url": f"{API_ROOT}/{quote(str(accession), safe='')}.json",
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
    parser.add_argument("--accession", default=DEFAULT_ACCESSION, help="UniProt accession such as P04637.")
    parser.add_argument("--feature-limit", type=int, default=8, help="Maximum number of representative features.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    payload = fetch_entry(args.accession)
    summary = build_summary(payload, feature_limit=max(1, args.feature_limit))
    write_json(summary, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
