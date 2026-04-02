#!/usr/bin/env python3
"""Generate or inspect a small FITS image with Astropy."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from astropy.io import fits


def build_toy_hdu() -> fits.PrimaryHDU:
    data = np.arange(25, dtype=np.float32).reshape((5, 5))
    hdu = fits.PrimaryHDU(data)
    hdu.header["FILTER"] = "r"
    hdu.header["EXPTIME"] = 30.0
    hdu.header["OBJECT"] = "ToyField"
    return hdu


def summarize_hdu(hdu: fits.PrimaryHDU, mode: str, source_path: Path | None) -> dict[str, object]:
    data = np.asarray(hdu.data, dtype=float)
    return {
        "mode": mode,
        "source_path": str(source_path) if source_path is not None else None,
        "shape": list(data.shape),
        "dtype": str(hdu.data.dtype),
        "mean_signal": round(float(data.mean()), 6),
        "std_signal": round(float(data.std()), 6),
        "min_signal": round(float(data.min()), 6),
        "max_signal": round(float(data.max()), 6),
        "filter": hdu.header.get("FILTER"),
        "exptime_s": round(float(hdu.header.get("EXPTIME")), 6),
        "object": hdu.header.get("OBJECT"),
        "bitpix": int(hdu.header.get("BITPIX")),
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
    parser.add_argument("--input", type=Path, default=None, help="Optional existing FITS file to summarize.")
    parser.add_argument("--fits-out", type=Path, default=None, help="Optional path to write the generated toy FITS.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    if args.input is not None and not args.input.exists():
        raise SystemExit(f"Input file not found: {args.input}")

    if args.input is not None:
        with fits.open(args.input) as hdul:
            payload = summarize_hdu(hdu=hdul[0], mode="input", source_path=args.input)
    else:
        hdu = build_toy_hdu()
        if args.fits_out is not None:
            args.fits_out.parent.mkdir(parents=True, exist_ok=True)
            hdu.writeto(args.fits_out, overwrite=True)
        payload = summarize_hdu(hdu=hdu, mode="generated_toy", source_path=args.fits_out)

    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
