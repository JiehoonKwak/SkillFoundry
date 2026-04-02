#!/usr/bin/env python3
"""Submit a short identifier list to Reactome Analysis Service and summarize top pathways."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


API_ROOT = "https://reactome.org/AnalysisService/identifiers/projection"
DEFAULT_IDENTIFIERS = ["BRCA1", "TP53"]
USER_AGENT = "SciSkillUniverse/0.2"


def normalize_identifiers(raw: str | None = None, input_path: Path | None = None) -> list[str]:
    parts: list[str] = []
    if raw:
        parts.extend(raw.replace("\n", ",").split(","))
    if input_path is not None:
        parts.extend(input_path.read_text(encoding="utf-8").replace("\n", ",").split(","))

    identifiers: list[str] = []
    seen: set[str] = set()
    for part in parts:
        token = part.strip()
        if not token:
            continue
        if token not in seen:
            identifiers.append(token)
            seen.add(token)
    if not identifiers:
        raise ValueError("At least one identifier is required.")
    return identifiers


def build_request_url(page_size: int, page: int) -> str:
    query = urlencode({"pageSize": page_size, "page": page})
    return f"{API_ROOT}?{query}"


def fetch_enrichment(identifiers: list[str], page_size: int, page: int) -> dict[str, Any]:
    payload = "\n".join(identifiers).encode("utf-8")
    request = Request(
        build_request_url(page_size=page_size, page=page),
        data=payload,
        headers={
            "Accept": "application/json",
            "Content-Type": "text/plain",
            "User-Agent": USER_AGENT,
        },
        method="POST",
    )
    try:
        with urlopen(request, timeout=30) as response:
            return json.load(response)
    except HTTPError as exc:
        message = exc.read().decode("utf-8", errors="replace").strip()
        raise SystemExit(f"Reactome enrichment request failed: HTTP {exc.code} {message}") from exc
    except URLError as exc:
        raise SystemExit(f"Reactome enrichment request failed: {exc.reason}") from exc


def _resource_summary(records: Any) -> list[dict[str, Any]]:
    if not isinstance(records, list):
        return []
    summaries: list[dict[str, Any]] = []
    for record in records:
        if not isinstance(record, dict):
            continue
        summaries.append(
            {
                "resource": record.get("resource"),
                "pathways": record.get("pathways"),
                "filtered": record.get("filtered"),
            }
        )
    return summaries


def _species_summary(records: Any) -> list[dict[str, Any]]:
    if not isinstance(records, list):
        return []
    summaries: list[dict[str, Any]] = []
    for record in records:
        if not isinstance(record, dict):
            continue
        summaries.append(
            {
                "db_id": record.get("dbId"),
                "tax_id": record.get("taxId"),
                "name": record.get("name"),
                "pathways": record.get("pathways"),
                "filtered": record.get("filtered"),
            }
        )
    return summaries


def _pathway_summary(record: dict[str, Any]) -> dict[str, Any]:
    entities = record.get("entities") or {}
    reactions = record.get("reactions") or {}
    species = record.get("species") or {}
    return {
        "stable_id": record.get("stId"),
        "db_id": record.get("dbId"),
        "name": record.get("name"),
        "species": species.get("name"),
        "species_tax_id": species.get("taxId"),
        "is_disease": bool(record.get("inDisease")),
        "is_llp": bool(record.get("llp")),
        "entities_found": entities.get("found"),
        "entities_total": entities.get("total"),
        "entities_ratio": entities.get("ratio"),
        "entities_p_value": entities.get("pValue"),
        "entities_fdr": entities.get("fdr"),
        "reactions_found": reactions.get("found"),
        "reactions_total": reactions.get("total"),
        "reactions_ratio": reactions.get("ratio"),
    }


def build_summary(payload: dict[str, Any], identifiers: list[str], page_size: int, page: int) -> dict[str, Any]:
    pathways = payload.get("pathways") or []
    identifiers_not_found = payload.get("identifiersNotFound")
    return {
        "identifiers": identifiers,
        "page_size": page_size,
        "page": page,
        "pathways_found": payload.get("pathwaysFound"),
        "identifiers_not_found": identifiers_not_found,
        "summary": payload.get("summary"),
        "resource_summary": _resource_summary(payload.get("resourceSummary")),
        "species_summary": _species_summary(payload.get("speciesSummary")),
        "warnings": payload.get("warnings") or [],
        "pathways": [_pathway_summary(record) for record in pathways if isinstance(record, dict)],
        "source_url": build_request_url(page_size=page_size, page=page),
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
    parser.add_argument(
        "--identifiers",
        default=",".join(DEFAULT_IDENTIFIERS),
        help="Comma-separated identifiers to analyze",
    )
    parser.add_argument("--input-file", type=Path, default=None, help="Optional text file of identifiers")
    parser.add_argument("--page-size", type=int, default=5, help="Number of pathways to return")
    parser.add_argument("--page", type=int, default=1, help="Page number for paginated results")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path")
    args = parser.parse_args()

    if args.page_size < 1:
        raise SystemExit("--page-size must be positive.")
    if args.page < 1:
        raise SystemExit("--page must be positive.")

    identifiers = normalize_identifiers(raw=args.identifiers, input_path=args.input_file)
    payload = fetch_enrichment(identifiers=identifiers, page_size=args.page_size, page=args.page)
    summary = build_summary(payload, identifiers=identifiers, page_size=args.page_size, page=args.page)
    write_json(summary, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
