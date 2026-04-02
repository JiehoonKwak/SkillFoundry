#!/usr/bin/env python3
"""Backwards-compatible wrapper for the bundled protein language model starter example."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


SKILL_ROOT = Path(__file__).resolve().parents[1]
PIPELINE = SKILL_ROOT / "scripts" / "run_protein_language_model_function_analysis.py"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=SKILL_ROOT / "assets" / "toy_protein_lm_summary.json")
    args = parser.parse_args()
    embeddings_out = args.out.with_suffix(".tsv")
    embeddings_out.parent.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(
        [
            sys.executable,
            str(PIPELINE),
            "--input",
            str(SKILL_ROOT / "examples" / "toy_sequences.fasta"),
            "--labels",
            str(SKILL_ROOT / "examples" / "toy_labels.tsv"),
            "--config",
            str(SKILL_ROOT / "examples" / "analysis_config.json"),
            "--embeddings-out",
            str(embeddings_out),
            "--summary-out",
            str(args.out),
        ],
        check=False,
    )
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
