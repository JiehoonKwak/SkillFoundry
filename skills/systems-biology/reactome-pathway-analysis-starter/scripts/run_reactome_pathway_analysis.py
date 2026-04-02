#!/usr/bin/env python3
"""Run a compact ranked pathway analysis summary against the Reactome Analysis Service."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


API_ROOT = "https://reactome.org/AnalysisService/identifiers/projection"
DEFAULT_IDENTIFIERS = ["BRCA1", "TP53", "EGFR"]
DEFAULT_FDR = 0.05
USER_AGENT = "SciSkillUniverse/0.2"
RETRYABLE_STATUS_CODES = {429, 500, 502, 503, 504}
ASSET_FALLBACK = Path(__file__).resolve().parents[1] / "assets" / "canonical_pathway_analysis.json"


class ReactomePathwayAnalysisError(RuntimeError):
    """Raised when the Reactome pathway analysis starter cannot obtain a usable result."""


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
        if not token or token in seen:
            continue
        identifiers.append(token)
        seen.add(token)
    if not identifiers:
        raise ValueError("At least one identifier is required.")
    return identifiers


def build_request_url(page_size: int, page: int) -> str:
    return f"{API_ROOT}?{urlencode({'pageSize': page_size, 'page': page})}"


def fetch_analysis(identifiers: list[str], page_size: int) -> dict[str, Any]:
    request = Request(
        build_request_url(page_size=page_size, page=1),
        data="\n".join(identifiers).encode("utf-8"),
        headers={
            "Accept": "application/json",
            "Content-Type": "text/plain",
            "User-Agent": USER_AGENT,
        },
        method="POST",
    )
    for attempt in range(3):
        try:
            with urlopen(request, timeout=30) as response:
                return json.load(response)
        except HTTPError as exc:
            message = exc.read().decode("utf-8", errors="replace").strip()
            if exc.code in RETRYABLE_STATUS_CODES and attempt < 2:
                time.sleep(attempt + 1)
                continue
            raise ReactomePathwayAnalysisError(
                f"Reactome pathway analysis failed: HTTP {exc.code} {message}"
            ) from exc
        except URLError as exc:
            if attempt < 2:
                time.sleep(attempt + 1)
                continue
            raise ReactomePathwayAnalysisError(
                f"Reactome pathway analysis failed: {exc.reason}"
            ) from exc
    raise ReactomePathwayAnalysisError("Reactome pathway analysis failed after retries.")


def _pathway_summary(pathway: dict[str, Any]) -> dict[str, Any]:
    entities = pathway.get("entities") or {}
    reactions = pathway.get("reactions") or {}
    species = pathway.get("species") or {}
    return {
        "stable_id": pathway.get("stId"),
        "name": pathway.get("name"),
        "species": species.get("name"),
        "entities_found": entities.get("found"),
        "entities_total": entities.get("total"),
        "entities_fdr": entities.get("fdr"),
        "entities_p_value": entities.get("pValue"),
        "reactions_found": reactions.get("found"),
        "reactions_total": reactions.get("total"),
    }


def build_summary(payload: dict[str, Any], identifiers: list[str], top_n: int, fdr_threshold: float) -> dict[str, Any]:
    pathways = [_pathway_summary(item) for item in (payload.get("pathways") or []) if isinstance(item, dict)]
    top_pathways = pathways[:top_n]
    fdr_values = [item["entities_fdr"] for item in pathways if item.get("entities_fdr") is not None]
    significant = [
        item
        for item in pathways
        if item.get("entities_fdr") is not None and float(item["entities_fdr"]) <= fdr_threshold
    ]
    return {
        "identifiers": identifiers,
        "pathways_found": payload.get("pathwaysFound"),
        "identifier_not_found_count": len(payload.get("identifiersNotFound") or []),
        "top_n": top_n,
        "fdr_threshold": fdr_threshold,
        "top_pathways": top_pathways,
        "significant_pathway_count": len(significant),
        "best_fdr": min(fdr_values) if fdr_values else None,
        "source_url": build_request_url(page_size=max(top_n, 10), page=1),
        "result_origin": "live_api",
    }


def load_asset_fallback(identifiers: list[str], reason: str) -> dict[str, Any] | None:
    if identifiers != DEFAULT_IDENTIFIERS or not ASSET_FALLBACK.exists():
        return None
    payload = json.loads(ASSET_FALLBACK.read_text(encoding="utf-8"))
    payload["result_origin"] = "asset_fallback"
    payload["fallback_reason"] = reason
    return payload


def write_json(payload: dict[str, Any], out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--identifiers", default=",".join(DEFAULT_IDENTIFIERS), help="Comma-separated identifiers.")
    parser.add_argument("--input-file", type=Path, default=None, help="Optional text file with identifiers.")
    parser.add_argument("--top-n", type=int, default=5, help="How many ranked pathways to keep.")
    parser.add_argument("--fdr-threshold", type=float, default=DEFAULT_FDR, help="Threshold for counting significant pathways.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    if args.top_n < 1:
        raise SystemExit("--top-n must be positive.")
    if args.fdr_threshold <= 0:
        raise SystemExit("--fdr-threshold must be positive.")

    identifiers = normalize_identifiers(raw=args.identifiers, input_path=args.input_file)
    try:
        payload = fetch_analysis(identifiers, page_size=max(args.top_n, 10))
        summary = build_summary(payload, identifiers=identifiers, top_n=args.top_n, fdr_threshold=args.fdr_threshold)
    except (ValueError, ReactomePathwayAnalysisError) as exc:
        fallback = load_asset_fallback(identifiers, str(exc))
        if fallback is None:
            raise SystemExit(str(exc)) from exc
        summary = fallback
    write_json(summary, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
