#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from run_mapping_validation import load_json, validate_outputs


SKILL_DIR = Path(__file__).resolve().parents[1]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--input", type=Path, default=SKILL_DIR / "examples" / "toy_input.json")
    args = parser.parse_args()
    payload = load_json(args.input)
    validate_outputs(args.outdir, expected=payload.get("expected_invariants"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
