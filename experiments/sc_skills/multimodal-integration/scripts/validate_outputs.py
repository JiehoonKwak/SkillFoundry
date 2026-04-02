#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from multimodal_integration_starter import validate_outputs


parser = argparse.ArgumentParser()
parser.add_argument("--outdir", type=Path, required=True)
args = parser.parse_args()
validate_outputs(Path(__file__).resolve().parents[1], args.outdir.resolve())
