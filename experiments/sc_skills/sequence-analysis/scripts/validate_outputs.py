#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from run_sequence_analysis import DEFAULT_INPUT, validate_outputs


parser = argparse.ArgumentParser()
parser.add_argument("--outdir", type=Path, required=True)
parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
args = parser.parse_args()
result = validate_outputs(args.outdir, input_path=args.input)
print(result)
