#!/usr/bin/env python3
"""Backwards-compatible wrapper for the promoted GWAS starter.

The benchmark suite still expects the old frontier-starter summary fields. Keep
emitting those fields while also exposing the richer QC summary from the
promoted starter implementation.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


SKILL_ROOT = Path(__file__).resolve().parents[1]
QC_SCRIPT = SKILL_ROOT / "scripts" / "run_gwas_summary_qc.py"
RESOURCE_CONTEXT = SKILL_ROOT / "examples" / "resource_context.json"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=SKILL_ROOT / "assets" / "toy_gwas_qc_summary.json")
    args = parser.parse_args()
    out_tsv = args.out.with_suffix(".tsv")
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(
        [
            sys.executable,
            str(QC_SCRIPT),
            "--input",
            str(SKILL_ROOT / "examples" / "toy_sumstats.tsv"),
            "--config",
            str(SKILL_ROOT / "examples" / "qc_config.json"),
            "--out-tsv",
            str(out_tsv),
            "--summary-out",
            str(args.out),
        ],
        check=False,
    )
    if completed.returncode != 0 or not args.out.exists():
        return completed.returncode

    payload = json.loads(args.out.read_text(encoding="utf-8"))
    context = json.loads(RESOURCE_CONTEXT.read_text(encoding="utf-8"))
    payload.update(
        {
            "skill_slug": context["skill_slug"],
            "domain_slug": context["domain_slug"],
            "leaf_slug": context["leaf_slug"],
            # Keep the legacy benchmark contract stable by surfacing the
            # canonical primary anchor first.
            "source_resource_ids": ["plink2-docs"],
            "starter_steps": list(context["starter_objectives"]),
            "promotion_checklist": [
                "Add a runnable example invocation to the skill entry point.",
                "Add or keep a repository-level smoke or regression test.",
                "Promote only after sandbox verification remains green.",
            ],
        }
    )
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
