#!/usr/bin/env python3
"""Search EBI BioSamples and return compact sample metadata summaries."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


API_ROOT = "https://www.ebi.ac.uk/biosamples/samples"
DEFAULT_TEXT = "breast cancer"
USER_AGENT = "SciSkillUniverse/0.2"


def normalize_text(text: str) -> str:
    normalized = " ".join(text.split())
    if not normalized:
        raise ValueError("Search text must not be empty.")
    return normalized


def normalize_size(size: int) -> int:
    if size <= 0:
        raise ValueError("Size must be a positive integer.")
    return size


def normalize_page(page: int) -> int:
    if page < 0:
        raise ValueError("Page must be zero or greater.")
    return page


def build_query_url(text: str, page: int, size: int) -> str:
    query = urlencode({"text": text, "page": page, "size": size})
    return f"{API_ROOT}?{query}"


def fetch_samples(text: str, page: int, size: int) -> dict[str, Any]:
    normalized_text = normalize_text(text)
    normalized_page = normalize_page(page)
    normalized_size = normalize_size(size)
    request = Request(
        build_query_url(normalized_text, normalized_page, normalized_size),
        headers={
            "Accept": "application/hal+json, application/json",
            "User-Agent": USER_AGENT,
        },
    )
    try:
        with urlopen(request, timeout=30) as response:
            payload = json.load(response)
    except HTTPError as exc:
        message = exc.read().decode("utf-8", errors="replace").strip()
        raise SystemExit(f"BioSamples request failed for text {text!r}: HTTP {exc.code} {message}") from exc
    except URLError as exc:
        raise SystemExit(f"BioSamples request failed for text {text!r}: {exc.reason}") from exc
    if not isinstance(payload, dict):
        raise SystemExit(f"Unexpected BioSamples response shape: expected dict, got {type(payload).__name__}")
    return payload


def _clean_text(value: Any) -> str | None:
    if isinstance(value, str):
        normalized = " ".join(value.split())
        return normalized or None
    return None


def _characteristic_values(characteristics: dict[str, Any], key: str, *, limit: int = 3) -> list[str]:
    raw_values = characteristics.get(key)
    if not isinstance(raw_values, list):
        return []
    values: list[str] = []
    for item in raw_values:
        text_value = None
        if isinstance(item, dict):
            text_value = _clean_text(item.get("text"))
        else:
            text_value = _clean_text(item)
        if text_value and text_value not in values:
            values.append(text_value)
        if len(values) >= limit:
            break
    return values


def summarize_sample(sample: dict[str, Any]) -> dict[str, Any]:
    characteristics = sample.get("characteristics")
    if not isinstance(characteristics, dict):
        characteristics = {}
    accession = _clean_text(sample.get("accession"))
    self_link = None
    links = sample.get("_links")
    if isinstance(links, dict) and isinstance(links.get("self"), dict):
        self_link = _clean_text(links["self"].get("href"))
    organism_values = _characteristic_values(characteristics, "organism", limit=3)
    disease_values = _characteristic_values(characteristics, "disease state", limit=3)
    organism_part_values = _characteristic_values(characteristics, "organism part", limit=3)
    return {
        "accession": accession,
        "name": _clean_text(sample.get("name")),
        "tax_id": sample.get("taxId"),
        "status": _clean_text(sample.get("status")),
        "release": _clean_text(sample.get("release")),
        "update": _clean_text(sample.get("update")),
        "organism": organism_values[0] if organism_values else None,
        "disease_state": disease_values[0] if disease_values else None,
        "organism_part": organism_part_values[0] if organism_part_values else None,
        "characteristics_keys": sorted(characteristics.keys())[:10],
        "sample_url": self_link or (f"https://www.ebi.ac.uk/biosamples/samples/{accession}" if accession else None),
    }


def build_summary(payload: dict[str, Any], text: str, page: int, size: int) -> dict[str, Any]:
    embedded = payload.get("_embedded")
    samples = embedded.get("samples") if isinstance(embedded, dict) else None
    if not isinstance(samples, list):
        raise SystemExit("Unexpected BioSamples response shape: missing _embedded.samples list.")
    page_info = payload.get("page")
    if not isinstance(page_info, dict):
        page_info = {}
    summaries = [summarize_sample(sample) for sample in samples]
    return {
        "text": normalize_text(text),
        "page": normalize_page(page),
        "size": normalize_size(size),
        "total_elements": page_info.get("totalElements"),
        "total_pages": page_info.get("totalPages"),
        "sample_count": len(summaries),
        "samples": summaries,
        "source_url": API_ROOT,
        "query_url": build_query_url(text, page, size),
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
    parser.add_argument("--text", default=DEFAULT_TEXT, help="Free-text BioSamples query.")
    parser.add_argument("--page", type=int, default=0, help="Zero-based results page.")
    parser.add_argument("--size", type=int, default=3, help="Maximum number of samples to summarize.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    try:
        payload = fetch_samples(args.text, args.page, args.size)
        summary = build_summary(payload, args.text, args.page, args.size)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    write_json(summary, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
