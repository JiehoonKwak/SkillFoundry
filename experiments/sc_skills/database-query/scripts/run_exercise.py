#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "scripts" / "run_database_query.py"
VALIDATOR = SKILL_DIR / "scripts" / "validate_outputs.py"
INPUT = SKILL_DIR / "examples" / "toy_input.json"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run the database-query toy exercise.")
    parser.add_argument("--outdir", type=Path, default=SKILL_DIR / "scratch" / "toy_run")
    args = parser.parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "--input", str(INPUT), "--outdir", str(args.outdir)],
        cwd=SKILL_DIR,
        check=False,
        text=True,
    )
    if completed.returncode != 0:
        return completed.returncode

    validated = subprocess.run(
        [sys.executable, str(VALIDATOR), "--input", str(INPUT), "--outdir", str(args.outdir)],
        cwd=SKILL_DIR,
        check=False,
        text=True,
    )
    return validated.returncode


if __name__ == "__main__":
    raise SystemExit(main())
