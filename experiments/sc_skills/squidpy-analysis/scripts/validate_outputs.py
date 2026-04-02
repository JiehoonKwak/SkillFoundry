#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from run_squidpy_analysis import SKILL_DIR, validate_outputs


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--input", type=Path, default=SKILL_DIR / "examples" / "toy_input.json")
    args = parser.parse_args()
    validate_outputs(SKILL_DIR, args.outdir, args.input)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
