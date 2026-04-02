#!/usr/bin/env python3
"""Detect near-duplicate resource records with datasketch."""

from __future__ import annotations

import argparse
import itertools
import json
import re
from pathlib import Path

from datasketch import MinHash, MinHashLSH


TOKEN_RE = re.compile(r"[a-z0-9]+")


def tokenize(text: str) -> set[str]:
    return set(TOKEN_RE.findall(text.lower()))


def load_resources(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def build_minhash(tokens: set[str], num_perm: int = 128) -> MinHash:
    mh = MinHash(num_perm=num_perm)
    for token in sorted(tokens):
        mh.update(token.encode("utf-8"))
    return mh


def run_dedup(input_path: Path, threshold: float) -> dict:
    resources = load_resources(input_path)
    token_sets = {
        item["resource_id"]: tokenize(f'{item["canonical_name"]} {item["summary"]}')
        for item in resources
    }
    minhashes = {resource_id: build_minhash(tokens) for resource_id, tokens in token_sets.items()}
    lsh = MinHashLSH(threshold=threshold, num_perm=128)
    for resource_id, mh in minhashes.items():
        lsh.insert(resource_id, mh)
    candidates = set()
    for resource_id, mh in minhashes.items():
        for hit in lsh.query(mh):
            if hit != resource_id:
                candidates.add(tuple(sorted((resource_id, hit))))
    pairs = []
    for left, right in sorted(candidates):
        score = minhashes[left].jaccard(minhashes[right])
        pairs.append({"left_resource_id": left, "right_resource_id": right, "jaccard_estimate": round(float(score), 6)})
    return {
        "input_path": str(input_path.resolve()),
        "resource_count": len(resources),
        "threshold": threshold,
        "candidate_pairs": pairs,
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
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    payload = run_dedup(args.input, args.threshold)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
