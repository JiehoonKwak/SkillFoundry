#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from run_panel_design import validate_outputs  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate panel-design starter outputs.")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    validate_outputs(Path(__file__).resolve().parents[1], args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
