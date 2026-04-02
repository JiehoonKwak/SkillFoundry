#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from run_spatial_domain_detection import SKILL_DIR, validate_outputs


parser = argparse.ArgumentParser()
parser.add_argument("--outdir", type=Path, required=True)
args = parser.parse_args()
validate_outputs(SKILL_DIR, args.outdir)
