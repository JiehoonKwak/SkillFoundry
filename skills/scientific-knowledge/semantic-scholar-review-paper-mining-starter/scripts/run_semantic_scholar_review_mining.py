#!/usr/bin/env python3
"""Identify likely review papers from local paper metadata."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "examples" / "paper_metadata.json"
REVIEW_PATTERNS = [
    re.compile(pattern, re.IGNORECASE)
    for pattern in (
        r"\breview\b",
        r"\bsurvey\b",
        r"\bmeta-analysis\b",
        r"\bsystematic review\b",
        r"\boverview\b",
    )
]


def review_hits(text: str) -> list[str]:
    return sorted({pattern.pattern for pattern in REVIEW_PATTERNS if pattern.search(text)})


def build_summary(input_path: Path, limit: int) -> dict[str, object]:
    papers = json.loads(input_path.read_text(encoding="utf-8"))
    if not isinstance(papers, list):
        raise SystemExit("Input JSON must contain a list of paper objects.")
    matches = []
    for paper in papers:
        text = " ".join(str(paper.get(field, "")) for field in ("title", "abstract", "publication_types"))
        hits = review_hits(text)
        if not hits:
            continue
        matches.append(
            {
                "paper_id": paper.get("paper_id"),
                "title": paper.get("title"),
                "year": paper.get("year"),
                "venue": paper.get("venue"),
                "citation_count": paper.get("citation_count", 0),
                "review_signals": hits,
                "review_score": len(hits) * 2 + min(int(paper.get("citation_count", 0)), 500) / 100,
            }
        )
    matches.sort(key=lambda item: (-float(item["review_score"]), -int(item["citation_count"]), -int(item["year"])))
    return {
        "input_path": str(input_path),
        "candidate_count": len(papers),
        "review_paper_count": len(matches),
        "review_papers": matches[:limit],
    }


def write_json(payload: dict[str, object], out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--limit", type=int, default=5)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()

    if not args.input.exists():
        raise SystemExit(f"Input JSON not found: {args.input}")
    if args.limit < 1:
        raise SystemExit("--limit must be >= 1")

    payload = build_summary(args.input, args.limit)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
