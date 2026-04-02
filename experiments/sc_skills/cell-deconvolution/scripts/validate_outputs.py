#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from run_cell_deconvolution import SKILL_DIR, validate_outputs


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate Tangram-style cell deconvolution starter outputs.")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    validate_outputs(SKILL_DIR, args.outdir.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
