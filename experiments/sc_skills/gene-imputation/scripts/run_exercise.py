#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / 'scripts' / 'run_gene_imputation.py'
VALIDATOR = SKILL_DIR / 'scripts' / 'validate_outputs.py'
INPUT = SKILL_DIR / 'examples' / 'toy_input.json'

parser = argparse.ArgumentParser()
parser.add_argument('--outdir', type=Path, default=SKILL_DIR / 'scratch' / 'toy_run')
args = parser.parse_args()
args.outdir.mkdir(parents=True, exist_ok=True)
completed = subprocess.run(
    [sys.executable, str(SCRIPT), '--input', str(INPUT), '--outdir', str(args.outdir)],
    cwd=SKILL_DIR.parents[2],
    check=False,
    text=True,
)
if completed.returncode != 0:
    raise SystemExit(completed.returncode)
validated = subprocess.run(
    [sys.executable, str(VALIDATOR), '--outdir', str(args.outdir)],
    cwd=SKILL_DIR.parents[2],
    check=False,
    text=True,
)
raise SystemExit(validated.returncode)
