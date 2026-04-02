#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from run_tissue_niche_annotation import validate_outputs


parser = argparse.ArgumentParser()
parser.add_argument('--outdir', type=Path, required=True)
args = parser.parse_args()
validate_outputs(Path(__file__).resolve().parents[1], args.outdir)
