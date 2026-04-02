#!/usr/bin/env python3
"""Compatibility wrapper for the canonical EBI Proteins summary script."""

from __future__ import annotations

import runpy
from pathlib import Path


if __name__ == "__main__":
    runpy.run_path(str(Path(__file__).with_name("fetch_protein_entry_summary.py")), run_name="__main__")
