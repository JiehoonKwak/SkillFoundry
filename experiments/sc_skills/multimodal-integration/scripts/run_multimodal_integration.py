#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from multimodal_integration_starter import run_cli


if __name__ == "__main__":
    raise SystemExit(run_cli(Path(__file__).resolve().parents[1]))
