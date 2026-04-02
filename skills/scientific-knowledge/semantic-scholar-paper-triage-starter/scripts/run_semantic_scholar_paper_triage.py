#!/usr/bin/env python3
"""Rank candidate papers for triage using deterministic heuristics."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "examples" / "candidate_papers.json"


def tokenize(text: str) -> set[str]:
    return {token for token in re.findall(r"[A-Za-z0-9]+", text.lower()) if len(token) > 2}


def score_paper(paper: dict[str, object], query_tokens: set[str]) -> tuple[float, dict[str, float]]:
    text = " ".join(
        str(paper.get(key, ""))
        for key in ("title", "abstract", "venue")
    )
    text_tokens = tokenize(text)
    overlap = len(query_tokens & text_tokens)
    citation_count = int(paper.get("citation_count", 0))
    year = int(paper.get("year", 0))
    overlap_score = overlap * 3.0
    citation_score = min(math.log1p(max(citation_count, 0)), 6.0)
    recency_score = max(min(year - 2018, 8), 0) * 0.2
    total = overlap_score + citation_score + recency_score
    return total, {
        "overlap_score": round(overlap_score, 3),
        "citation_score": round(citation_score, 3),
        "recency_score": round(recency_score, 3),
    }


def build_summary(input_path: Path, query: str, limit: int) -> dict[str, object]:
    papers = json.loads(input_path.read_text(encoding="utf-8"))
    if not isinstance(papers, list) or not papers:
        raise SystemExit("Input JSON must contain a non-empty list of paper objects.")
    query_tokens = tokenize(query)
    ranked = []
    for paper in papers:
        score, components = score_paper(paper, query_tokens)
        ranked.append(
            {
                "paper_id": paper.get("paper_id"),
                "title": paper.get("title"),
                "year": paper.get("year"),
                "venue": paper.get("venue"),
                "citation_count": paper.get("citation_count"),
                "triage_score": round(score, 3),
                "score_components": components,
            }
        )
    ranked.sort(key=lambda item: (-float(item["triage_score"]), -int(item["citation_count"]), -int(item["year"])))
    return {
        "query": query,
        "input_path": str(input_path),
        "candidate_count": len(papers),
        "top_candidates": ranked[:limit],
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
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Input JSON file with candidate paper metadata.")
    parser.add_argument("--query", required=True, help="Free-text triage query used for lexical overlap scoring.")
    parser.add_argument("--limit", type=int, default=5, help="Number of ranked papers to keep.")
    parser.add_argument("--out", type=Path, default=None, help="Optional output JSON path.")
    args = parser.parse_args()

    if not args.input.exists():
        raise SystemExit(f"Input JSON not found: {args.input}")
    if args.limit < 1:
        raise SystemExit("--limit must be >= 1")

    payload = build_summary(args.input, args.query, args.limit)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
