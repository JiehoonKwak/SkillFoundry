#!/usr/bin/env python3
"""Generate or inspect a small NIfTI image with NiBabel."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import nibabel as nib
import numpy as np


def build_toy_image() -> nib.Nifti1Image:
    array = np.arange(4 * 4 * 3 * 2, dtype=np.float32).reshape((4, 4, 3, 2)) / 10.0
    affine = np.array(
        [
            [2.0, 0.0, 0.0, -4.0],
            [0.0, 2.0, 0.0, -4.0],
            [0.0, 0.0, 2.0, -3.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        dtype=float,
    )
    image = nib.Nifti1Image(array, affine)
    image.header.set_xyzt_units("mm", "sec")
    image.header["pixdim"][4] = 1.5
    return image


def summarize_image(image: nib.Nifti1Image, mode: str, source_path: Path | None) -> dict[str, object]:
    data = image.get_fdata()
    spatial_units, time_units = image.header.get_xyzt_units()
    return {
        "mode": mode,
        "source_path": str(source_path) if source_path is not None else None,
        "shape": list(image.shape),
        "zooms": [round(float(value), 6) for value in image.header.get_zooms()],
        "dtype": str(image.get_data_dtype()),
        "mean_intensity": round(float(data.mean()), 6),
        "first_voxel": round(float(data.flat[0]), 6),
        "last_voxel": round(float(data.flat[-1]), 6),
        "affine_determinant": round(float(np.linalg.det(image.affine[:3, :3])), 6),
        "spatial_units": spatial_units,
        "time_units": time_units,
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
    parser.add_argument("--input", type=Path, default=None, help="Optional existing NIfTI file to summarize.")
    parser.add_argument("--nifti-out", type=Path, default=None, help="Optional path to write the generated toy NIfTI.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    if args.input is not None and not args.input.exists():
        raise SystemExit(f"Input file not found: {args.input}")

    if args.input is not None:
        image = nib.load(args.input)
        payload = summarize_image(image=image, mode="input", source_path=args.input)
    else:
        image = build_toy_image()
        if args.nifti_out is not None:
            args.nifti_out.parent.mkdir(parents=True, exist_ok=True)
            nib.save(image, args.nifti_out)
        payload = summarize_image(image=image, mode="generated_toy", source_path=args.nifti_out)

    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
