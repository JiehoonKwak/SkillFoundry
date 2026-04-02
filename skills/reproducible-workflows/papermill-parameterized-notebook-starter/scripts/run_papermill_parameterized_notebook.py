#!/usr/bin/env python3
"""Execute a parameterized notebook with Papermill and summarize the output."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import nbformat
import papermill as pm


DEFAULT_INPUT = Path(__file__).resolve().parents[1] / "examples" / "toy_parameters.ipynb"


def extract_result(output_notebook: Path) -> tuple[dict, bool, int]:
    notebook = nbformat.read(output_notebook, as_version=4)
    injected_parameters_present = False
    executed_code_cells = 0
    result: dict | None = None

    for cell in notebook.cells:
        if cell.get("cell_type") == "code":
            if cell.get("execution_count") is not None:
                executed_code_cells += 1
            tags = cell.get("metadata", {}).get("tags", [])
            if "injected-parameters" in tags:
                injected_parameters_present = True
            for output in cell.get("outputs", []):
                if output.get("output_type") != "stream":
                    continue
                text = output.get("text", "").strip()
                if not text:
                    continue
                try:
                    maybe_result = json.loads(text)
                except json.JSONDecodeError:
                    continue
                if isinstance(maybe_result, dict) and {"x", "y", "sum", "product"} <= set(maybe_result):
                    result = maybe_result
    if result is None:
        raise ValueError("Could not extract the final JSON result from the executed notebook.")
    return result, injected_parameters_present, executed_code_cells


def run_notebook(input_path: Path, output_notebook: Path, x: int, y: int) -> dict:
    output_notebook.parent.mkdir(parents=True, exist_ok=True)
    pm.execute_notebook(
        str(input_path),
        str(output_notebook),
        parameters={"x": x, "y": y},
        kernel_name="python3",
        log_output=False,
    )
    result, injected_parameters_present, executed_code_cells = extract_result(output_notebook)
    return {
        "input_notebook": str(input_path),
        "output_notebook": str(output_notebook),
        "parameters": {"x": int(x), "y": int(y)},
        "result": {
            "x": int(result["x"]),
            "y": int(result["y"]),
            "sum": int(result["sum"]),
            "product": int(result["product"]),
        },
        "executed_code_cell_count": executed_code_cells,
        "injected_parameters_present": injected_parameters_present,
        "output_notebook_size_bytes": int(output_notebook.stat().st_size),
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
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Notebook path to execute.")
    parser.add_argument("--output-notebook", type=Path, required=True, help="Executed notebook output path.")
    parser.add_argument("--summary-out", type=Path, default=None, help="Optional JSON summary output path.")
    parser.add_argument("--x", type=int, default=5, help="First integer parameter.")
    parser.add_argument("--y", type=int, default=7, help="Second integer parameter.")
    args = parser.parse_args()

    if not args.input.exists():
        raise SystemExit(f"Notebook not found: {args.input}")

    payload = run_notebook(args.input, args.output_notebook, args.x, args.y)
    write_json(payload, args.summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
