#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve().with_name("run_cellphonedb_analysis.py")

spec = importlib.util.spec_from_file_location("cellphonedb_starter", SCRIPT_PATH)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)

if __name__ == "__main__":
    raise SystemExit(module.main())
