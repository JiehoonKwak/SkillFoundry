#!/usr/bin/env python3
"""Generate a toy raster and summarize a windowed Rasterio preprocessing pass."""

from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path


def load_stack():
    try:
        import numpy as np
        import rasterio
        from rasterio.enums import Resampling
        from rasterio.transform import from_origin
        from rasterio.windows import Window
    except ImportError as exc:  # pragma: no cover - exercised via runtime tests
        raise SystemExit(
            "This script requires the geospatial prefix with rasterio installed at slurm/envs/geospatial."
        ) from exc
    return np, rasterio, Resampling, from_origin, Window


def write_toy_raster(path: Path) -> None:
    np, rasterio, _, from_origin, _ = load_stack()
    data = np.arange(1, 17, dtype="float32").reshape(4, 4)
    profile = {
        "driver": "GTiff",
        "height": 4,
        "width": 4,
        "count": 1,
        "dtype": "float32",
        "crs": "EPSG:4326",
        "transform": from_origin(0, 4, 1, 1),
    }
    with rasterio.open(path, "w", **profile) as dataset:
        dataset.write(data, 1)


def summarize_raster(path: Path) -> dict:
    np, rasterio, Resampling, _, Window = load_stack()
    with rasterio.open(path) as dataset:
        window = Window(1, 1, 2, 2)
        window_values = dataset.read(1, window=window)
        resampled = dataset.read(1, out_shape=(1, 2, 2), resampling=Resampling.average)
        return {
            "input_shape": [int(dataset.height), int(dataset.width)],
            "crs": str(dataset.crs),
            "bounds": {
                "left": float(dataset.bounds.left),
                "bottom": float(dataset.bounds.bottom),
                "right": float(dataset.bounds.right),
                "top": float(dataset.bounds.top),
            },
            "window": {
                "col_off": int(window.col_off),
                "row_off": int(window.row_off),
                "width": int(window.width),
                "height": int(window.height),
            },
            "window_values": [[round(float(value), 6) for value in row] for row in window_values.tolist()],
            "window_mean": round(float(np.mean(window_values)), 6),
            "resampled_shape": [int(value) for value in resampled.shape],
            "resampled_values": [[round(float(value), 6) for value in row] for row in resampled.tolist()],
            "resampled_mean": round(float(np.mean(resampled)), 6),
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
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON summary output path.")
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as tmp_dir:
        raster_path = Path(tmp_dir) / "toy_raster.tif"
        write_toy_raster(raster_path)
        payload = summarize_raster(raster_path)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
