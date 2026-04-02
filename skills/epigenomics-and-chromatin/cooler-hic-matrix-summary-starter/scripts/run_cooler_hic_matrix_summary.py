#!/usr/bin/env python3
"""Generate or summarize a tiny Hi-C matrix stored in Cooler format."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import cooler
import pandas as pd


ROOT = Path(__file__).resolve().parents[4]


def display_path(path: Path | None) -> str | None:
    if path is None:
        return None
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def build_toy_cooler(cooler_out: Path) -> None:
    bins = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1"],
            "start": [0, 1000, 2000],
            "end": [1000, 2000, 3000],
        }
    )
    pixels = pd.DataFrame(
        {
            "bin1_id": [0, 0, 1, 1, 2],
            "bin2_id": [0, 1, 1, 2, 2],
            "count": [10, 4, 12, 5, 7],
        }
    )
    cooler_out.parent.mkdir(parents=True, exist_ok=True)
    if cooler_out.exists():
        cooler_out.unlink()
    cooler.create_cooler(str(cooler_out), bins=bins, pixels=pixels, ordered=True, dtypes={"count": "int32"})


def summarize_cooler(cool_path: Path) -> dict[str, object]:
    clr = cooler.Cooler(str(cool_path))
    bins = clr.bins()[:]
    pixels = clr.pixels()[:]
    return {
        "cooler_path": display_path(cool_path),
        "bin_count": int(len(bins)),
        "pixel_count": int(len(pixels)),
        "matrix_shape": [int(value) for value in clr.shape],
        "chromnames": [str(name) for name in clr.chromnames],
        "total_contact_count": int(pixels["count"].sum()),
        "first_pixels": pixels.head(3).to_dict(orient="records"),
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
    parser.add_argument("--input", type=Path, default=None, help="Optional existing .cool input.")
    parser.add_argument(
        "--cooler-out",
        type=Path,
        default=Path("scratch/epigenomics/toy_contacts.cool"),
        help="Output path for the generated toy .cool file when --input is omitted.",
    )
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON summary path.")
    args = parser.parse_args()

    if args.input is not None:
        cool_path = args.input if args.input.is_absolute() else ROOT / args.input
        if not cool_path.exists():
            raise SystemExit(f"Cooler input not found: {cool_path}")
    else:
        cool_path = args.cooler_out if args.cooler_out.is_absolute() else ROOT / args.cooler_out
        build_toy_cooler(cool_path.resolve())

    out_path = args.out if args.out is None or args.out.is_absolute() else ROOT / args.out
    payload = summarize_cooler(cool_path.resolve())
    write_json(payload, out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
