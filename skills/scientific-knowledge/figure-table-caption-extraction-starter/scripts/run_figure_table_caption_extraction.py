#!/usr/bin/env python3
"""Extract figure and table captions from local text."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "examples" / "paper_excerpt.txt"
CAPTION_PATTERN = re.compile(r"^(Figure|Fig\.|Table)\s+(\d+)[\.:]?\s*(.+)$", re.IGNORECASE)


def build_summary(input_path: Path) -> dict[str, object]:
    figure_captions = []
    table_captions = []
    for raw_line in input_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        match = CAPTION_PATTERN.match(line)
        if not match:
            continue
        label, number, caption = match.groups()
        record = {"label": label, "number": number, "caption": caption}
        if label.lower().startswith("table"):
            table_captions.append(record)
        else:
            figure_captions.append(record)
    return {
        "input_path": str(input_path),
        "figure_caption_count": len(figure_captions),
        "table_caption_count": len(table_captions),
        "figure_captions": figure_captions,
        "table_captions": table_captions,
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
