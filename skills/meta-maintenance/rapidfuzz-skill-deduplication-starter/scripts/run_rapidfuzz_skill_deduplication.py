#!/usr/bin/env python3
"""Detect near-duplicate skill names with RapidFuzz."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
from pathlib import Path

from rapidfuzz import fuzz


def load_skills(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run_dedup(input_path: Path, threshold: int) -> dict:
    skills = load_skills(input_path)
    candidates = []
    for left, right in itertools.combinations(skills, 2):
        score = int(round(fuzz.token_sort_ratio(left["name"], right["name"])))
        if score >= threshold:
            candidates.append(
                {
                    "left_slug": left["slug"],
                    "right_slug": right["slug"],
                    "score": score,
                }
            )
    candidates.sort(key=lambda item: (-item["score"], item["left_slug"], item["right_slug"]))
    return {
        "input_path": str(input_path.resolve()),
        "threshold": threshold,
        "skill_count": len(skills),
        "candidate_pairs": candidates,
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
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--threshold", type=int, default=85)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    payload = run_dedup(args.input, args.threshold)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
