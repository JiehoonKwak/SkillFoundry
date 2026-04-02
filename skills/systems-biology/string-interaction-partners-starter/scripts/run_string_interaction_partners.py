#!/usr/bin/env python3
"""Query STRING interaction partners for one identifier and summarize the result."""

from __future__ import annotations

import argparse
import json
import socket
import time
import urllib.parse
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


API_ROOT = "https://string-db.org/api/json/interaction_partners"
DEFAULT_IDENTIFIER = "TP53"
DEFAULT_SPECIES = "9606"
DEFAULT_LIMIT = 5
DEFAULT_REQUIRED_SCORE = 700
USER_AGENT = "SciSkillUniverse/0.2"
ASSET_FALLBACK = Path(__file__).resolve().parents[1] / "assets" / "tp53_interaction_partners.json"
RETRYABLE_STATUS_CODES = {408, 429, 500, 502, 503, 504}


class StringInteractionError(RuntimeError):
    """Raised when the STRING API cannot be queried reliably."""


def normalize_identifier(identifier: str) -> str:
    normalized = identifier.strip()
    if not normalized:
        raise ValueError("Identifier must not be empty.")
    return normalized


def normalize_species(species: str) -> str:
    normalized = species.strip()
    if not normalized:
        raise ValueError("Species must not be empty.")
    return normalized


def read_identifier(identifier: str | None, identifier_file: Path | None) -> str:
    if identifier_file is not None:
        if not identifier_file.exists():
            raise ValueError(f"Identifier file does not exist: {identifier_file}")
        lines = [line.strip() for line in identifier_file.read_text(encoding="utf-8").splitlines() if line.strip()]
        if not lines:
            raise ValueError("Identifier file must contain at least one non-empty line.")
        return normalize_identifier(lines[0])
    if identifier is None:
        return DEFAULT_IDENTIFIER
    return normalize_identifier(identifier)


def fetch_interaction_partners(identifier: str, species: str, limit: int, required_score: int) -> list[dict[str, Any]]:
    payload = urllib.parse.urlencode(
        {
            "identifiers": normalize_identifier(identifier),
            "species": normalize_species(species),
            "limit": str(limit),
            "required_score": str(required_score),
            "caller_identity": USER_AGENT,
        }
    ).encode("utf-8")
    request = Request(
        API_ROOT,
        data=payload,
        headers={
            "Accept": "application/json",
            "Content-Type": "application/x-www-form-urlencoded",
            "User-Agent": USER_AGENT,
        },
    )
    for attempt in range(4):
        try:
            with urlopen(request, timeout=30) as response:
                result = json.load(response)
            if not isinstance(result, list):
                raise StringInteractionError(f"Unexpected STRING payload type: {type(result).__name__}")
            return result
        except HTTPError as exc:
            message = exc.read().decode("utf-8", errors="replace").strip()
            if exc.code in RETRYABLE_STATUS_CODES and attempt < 3:
                delay = float(exc.headers.get("Retry-After") or (attempt + 1))
                time.sleep(delay)
                continue
            raise StringInteractionError(
                f"STRING interaction_partners request failed for {identifier}: HTTP {exc.code} {message}"
            ) from exc
        except URLError as exc:
            if attempt < 3:
                time.sleep(float(attempt + 1))
                continue
            raise StringInteractionError(f"STRING interaction_partners request failed for {identifier}: {exc.reason}") from exc
        except (TimeoutError, socket.timeout) as exc:
            if attempt < 3:
                time.sleep(float(attempt + 1))
                continue
            raise StringInteractionError(f"STRING interaction_partners request failed for {identifier}: {exc}") from exc


def build_summary(records: list[dict[str, Any]], identifier: str, species: str, limit: int, required_score: int) -> dict[str, Any]:
    if not records:
        raise StringInteractionError(f"STRING returned no interaction partners for {identifier} in species {species}.")
    query_string_id = records[0].get("stringId_A")
    query_preferred_name = records[0].get("preferredName_A") or identifier
    top_partners = []
    scores = []
    for record in records[:limit]:
        score = float(record.get("score", 0.0))
        scores.append(score)
        top_partners.append(
            {
                "partner_string_id": record.get("stringId_B"),
                "partner_preferred_name": record.get("preferredName_B"),
                "score": round(score, 6),
                "evidence_scores": {
                    "neighborhood": round(float(record.get("nscore", 0.0)), 6),
                    "fusion": round(float(record.get("fscore", 0.0)), 6),
                    "cooccurrence": round(float(record.get("pscore", 0.0)), 6),
                    "coexpression": round(float(record.get("ascore", 0.0)), 6),
                    "experiments": round(float(record.get("escore", 0.0)), 6),
                    "database": round(float(record.get("dscore", 0.0)), 6),
                    "textmining": round(float(record.get("tscore", 0.0)), 6),
                },
            }
        )
    return {
        "query_identifier": identifier,
        "query_preferred_name": query_preferred_name,
        "query_string_id": query_string_id,
        "species": normalize_species(species),
        "limit": limit,
        "required_score": required_score,
        "partner_count": len(top_partners),
        "top_partners": top_partners,
        "score_max": round(max(scores), 6),
        "score_min": round(min(scores), 6),
        "source_url": API_ROOT,
        "result_origin": "live_api",
    }


def load_asset_fallback(identifier: str, species: str, limit: int, required_score: int, reason: str) -> dict[str, Any] | None:
    if (
        normalize_identifier(identifier) != DEFAULT_IDENTIFIER
        or normalize_species(species) != DEFAULT_SPECIES
        or limit != DEFAULT_LIMIT
        or required_score != DEFAULT_REQUIRED_SCORE
        or not ASSET_FALLBACK.exists()
    ):
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
    parser.add_argument("--identifier", default=None, help="One STRING-supported identifier, for example TP53.")
    parser.add_argument("--identifier-file", type=Path, default=None, help="Optional one-line text file with the identifier.")
    parser.add_argument("--species", default=DEFAULT_SPECIES, help="NCBI species identifier, default 9606.")
    parser.add_argument("--limit", type=int, default=DEFAULT_LIMIT, help="Maximum number of partners to include.")
    parser.add_argument(
        "--required-score",
        type=int,
        default=DEFAULT_REQUIRED_SCORE,
        help="Minimum STRING score threshold as an integer between 0 and 1000.",
    )
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    try:
        identifier = read_identifier(args.identifier, args.identifier_file)
        species = normalize_species(args.species)
        if args.limit < 1:
            raise ValueError("--limit must be positive.")
        if args.required_score < 0 or args.required_score > 1000:
            raise ValueError("--required-score must be between 0 and 1000.")
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    try:
        records = fetch_interaction_partners(identifier, species, args.limit, args.required_score)
        payload = build_summary(records, identifier, species, args.limit, args.required_score)
    except StringInteractionError as exc:
        fallback = load_asset_fallback(
            identifier,
            species,
            args.limit,
            args.required_score,
            str(exc),
        )
        if fallback is None:
            raise SystemExit(str(exc)) from exc
        payload = fallback

    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
