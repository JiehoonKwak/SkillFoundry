#!/usr/bin/env python3
"""Extract dataset, code, and package links from local paper text."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "examples" / "paper_text.md"
URL_PATTERN = re.compile(r"https?://[^\s)>\"]+")


def classify_url(url: str) -> str:
    lowered = url.lower()
    if any(token in lowered for token in ("github.com", "gitlab.com", "bitbucket.org")):
        return "code"
    if any(token in lowered for token in ("zenodo.org", "figshare.com", "kaggle.com", "huggingface.co/datasets", "dryad", "ebi.ac.uk/biostudies", "dataset")):
        return "dataset"
    if any(token in lowered for token in ("pypi.org", "bioconductor.org", "cran.r-project.org", "conda-forge")):
        return "package"
    return "other"


def build_summary(input_path: Path) -> dict[str, object]:
    text = input_path.read_text(encoding="utf-8")
    urls = sorted(dict.fromkeys(URL_PATTERN.findall(text)))
    classified = {"dataset": [], "code": [], "package": [], "other": []}
    for url in urls:
        classified[classify_url(url)].append(url)
    return {
        "input_path": str(input_path),
        "url_count": len(urls),
        "dataset_links": classified["dataset"],
        "code_links": classified["code"],
        "package_links": classified["package"],
        "other_links": classified["other"],
    }


def write_json(payload: dict[str, object], out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    if not args.input.exists():
        raise SystemExit(f"Input text not found: {args.input}")
    write_json(build_summary(args.input), args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
