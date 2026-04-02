#!/usr/bin/env python3
"""Render a tiny Jupyter notebook to standalone HTML with Quarto."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
REPORTING_PREFIX = ROOT / "slurm" / "envs" / "reporting"
QUARTO_BIN = REPORTING_PREFIX / "bin" / "quarto"
DEFAULT_INPUT = Path(__file__).resolve().parents[1] / "examples" / "toy_report.ipynb"


def quarto_env() -> dict[str, str]:
    env = os.environ.copy()
    prefix = str(REPORTING_PREFIX)
    env["PATH"] = f"{REPORTING_PREFIX / 'bin'}:{env.get('PATH', '')}"
    env["DENO_INSTALL_ROOT"] = prefix
    env["DENO_DOM_PLUGIN"] = str(REPORTING_PREFIX / "lib" / "deno_dom.so")
    env["DENO_DOM_VERSION"] = "0.1.41-alpha"
    env["QUARTO_DENO"] = str(REPORTING_PREFIX / "bin" / "deno")
    env["QUARTO_DENO_DOM"] = str(REPORTING_PREFIX / "lib" / "deno_dom.so")
    env["QUARTO_PANDOC"] = str(REPORTING_PREFIX / "bin" / "pandoc")
    env["QUARTO_ESBUILD"] = str(REPORTING_PREFIX / "bin" / "esbuild")
    env["QUARTO_TYPST"] = str(REPORTING_PREFIX / "bin" / "typst")
    env["QUARTO_DART_SASS"] = str(REPORTING_PREFIX / "bin" / "sass")
    env["QUARTO_SHARE_PATH"] = str(REPORTING_PREFIX / "share" / "quarto")
    env["QUARTO_CONDA_PREFIX"] = prefix
    return env


def render_notebook(input_path: Path, html_out: Path) -> dict:
    if not QUARTO_BIN.exists():
        raise SystemExit("Quarto is not installed at slurm/envs/reporting/bin/quarto")

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_root = Path(tmp_dir)
        notebook_path = tmp_root / "report.ipynb"
        output_dir = tmp_root / "out"
        output_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(input_path, notebook_path)
        subprocess.run(
            [
                str(QUARTO_BIN),
                "render",
                notebook_path.name,
                "--to",
                "html",
                "--output",
                "report.html",
                "--output-dir",
                "out",
                "--execute",
                "--no-cache",
            ],
            cwd=tmp_root,
            check=True,
            capture_output=True,
            text=True,
            timeout=300,
            env=quarto_env(),
        )
        rendered_html = output_dir / "report.html"
        html_out.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(rendered_html, html_out)

    html_text = html_out.read_text(encoding="utf-8")
    return {
        "input_notebook": str(input_path),
        "html_output": str(html_out),
        "title_present": "Toy Report" in html_text,
        "executed_output_present": "\"total\": 20" in html_text and "\"mean\": 5.0" in html_text,
        "html_size_bytes": int(html_out.stat().st_size),
    }


def write_json(payload: dict, out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Notebook path to render.")
    parser.add_argument("--html-out", type=Path, required=True, help="HTML output path.")
    parser.add_argument("--summary-out", type=Path, default=None, help="Optional JSON summary output path.")
    args = parser.parse_args()

    if not args.input.exists():
        raise SystemExit(f"Notebook not found: {args.input}")

    payload = render_notebook(input_path=args.input, html_out=args.html_out)
    write_json(payload, args.summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
