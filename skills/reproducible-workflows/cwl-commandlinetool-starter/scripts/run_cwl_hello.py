#!/usr/bin/env python3
"""Run the minimal local CWL CommandLineTool example and summarize its outputs."""

from __future__ import annotations

import argparse
import json
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SKILL_DIR = Path(__file__).resolve().parents[1]
CWLTOOL = ROOT / "slurm" / "envs" / "workflow-languages" / "bin" / "cwltool"
EXAMPLE_CWL = SKILL_DIR / "examples" / "hello.cwl"


def build_job_file(workspace: Path, message: str) -> Path:
    job_path = workspace / "hello-job.yml"
    job_path.write_text(f"message: {json.dumps(message)}\n", encoding="utf-8")
    return job_path


def run_command(message: str, workspace: Path) -> dict:
    if not CWLTOOL.exists():
        raise SystemExit(f"Missing cwltool executable: {CWLTOOL}")
    workspace.mkdir(parents=True, exist_ok=True)
    out_dir = workspace / "outputs"
    summary_json = workspace / "cwltool-summary.json"
    job_path = build_job_file(workspace, message)

    validate = subprocess.run(
        [str(CWLTOOL), "--validate", str(EXAMPLE_CWL)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
        timeout=120,
    )
    subprocess.run(
        [
            str(CWLTOOL),
            "--quiet",
            "--no-container",
            "--outdir",
            str(out_dir),
            "--write-summary",
            str(summary_json),
            str(EXAMPLE_CWL),
            str(job_path),
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
        timeout=120,
    )
    summary = json.loads(summary_json.read_text(encoding="utf-8"))
    output_path = Path(summary["greeting_file"]["location"].removeprefix("file://"))
    return {
        "message": message,
        "cwltool_version": subprocess.run(
            [str(CWLTOOL), "--version"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        ).stdout.strip().split()[-1],
        "validated_document": str(EXAMPLE_CWL.relative_to(ROOT)),
        "validation_message": validate.stderr.strip() or validate.stdout.strip(),
        "output_file": str(output_path),
        "output_text": output_path.read_text(encoding="utf-8").strip(),
        "raw_summary": summary,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message", default="hello from cwl", help="Greeting text to write")
    parser.add_argument("--workspace", type=Path, default=None, help="Optional workspace directory")
    parser.add_argument("--summary-out", type=Path, default=None, help="Optional JSON summary output")
    args = parser.parse_args()

    if args.workspace is None:
        with tempfile.TemporaryDirectory(prefix="cwl-hello-") as tmp_dir:
            payload = run_command(args.message, Path(tmp_dir))
    else:
        payload = run_command(args.message, args.workspace)

    text = json.dumps(payload, indent=2, sort_keys=True)
    if args.summary_out is not None:
        args.summary_out.parent.mkdir(parents=True, exist_ok=True)
        args.summary_out.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
