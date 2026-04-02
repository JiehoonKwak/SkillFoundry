#!/usr/bin/env python3
"""Compute a deterministic sourmash signature comparison for two toy FASTA files."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


DEFAULT_QUERY = Path(__file__).resolve().parents[1] / "examples" / "query.fa"
DEFAULT_REFERENCE = Path(__file__).resolve().parents[1] / "examples" / "reference.fa"


def load_stack():
    try:
        from sourmash import MinHash
    except ImportError as exc:  # pragma: no cover - exercised via runtime tests
        raise SystemExit(
            "This script requires the metagenomics prefix at slurm/envs/metagenomics."
        ) from exc
    return MinHash


def load_first_sequence(path: Path) -> tuple[str, str]:
    header = None
    chunks: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                break
            header = line[1:].strip() or path.stem
            continue
        chunks.append(line.upper())
    if header is None or not chunks:
        raise SystemExit(f"No FASTA record found in {path}")
    return header, "".join(chunks)


def build_minhash(sequence: str, ksize: int, scaled: int):
    MinHash = load_stack()
    mh = MinHash(n=0, ksize=ksize, scaled=scaled)
    mh.add_sequence(sequence, force=True)
    return mh


def compare_sequences(query_name: str, query_sequence: str, reference_name: str, reference_sequence: str, ksize: int, scaled: int) -> dict:
    query_mh = build_minhash(query_sequence, ksize=ksize, scaled=scaled)
    reference_mh = build_minhash(reference_sequence, ksize=ksize, scaled=scaled)
    query_hashes = set(query_mh.hashes)
    reference_hashes = set(reference_mh.hashes)
    shared_hashes = len(query_hashes & reference_hashes)
    return {
        "query_name": query_name,
        "reference_name": reference_name,
        "ksize": int(ksize),
        "scaled": int(scaled),
        "query_length": len(query_sequence),
        "reference_length": len(reference_sequence),
        "query_hash_count": len(query_hashes),
        "reference_hash_count": len(reference_hashes),
        "shared_hash_count": shared_hashes,
        "jaccard_similarity": round(float(query_mh.jaccard(reference_mh)), 6),
        "query_containment_in_reference": round(float(query_mh.contained_by(reference_mh)), 6),
        "reference_containment_in_query": round(float(reference_mh.contained_by(query_mh)), 6),
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
    parser.add_argument("--query-fasta", type=Path, default=DEFAULT_QUERY, help="Query FASTA path.")
    parser.add_argument("--reference-fasta", type=Path, default=DEFAULT_REFERENCE, help="Reference FASTA path.")
    parser.add_argument("--ksize", type=int, default=7, help="sourmash k-mer size.")
    parser.add_argument("--scaled", type=int, default=1, help="Scaled MinHash factor.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    if args.ksize <= 0:
        raise SystemExit("ksize must be a positive integer.")
    if args.scaled <= 0:
        raise SystemExit("scaled must be a positive integer.")

    query_name, query_sequence = load_first_sequence(args.query_fasta)
    reference_name, reference_sequence = load_first_sequence(args.reference_fasta)
    payload = compare_sequences(
        query_name=query_name,
        query_sequence=query_sequence,
        reference_name=reference_name,
        reference_sequence=reference_sequence,
        ksize=args.ksize,
        scaled=args.scaled,
    )
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
