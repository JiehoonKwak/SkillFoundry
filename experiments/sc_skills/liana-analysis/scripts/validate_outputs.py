#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from run_liana_analysis import SKILL_DIR, validate_outputs


parser = argparse.ArgumentParser()
parser.add_argument("--outdir", type=Path, required=True)
args = parser.parse_args()
validate_outputs(SKILL_DIR, args.outdir)
