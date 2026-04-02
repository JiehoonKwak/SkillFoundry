#!/usr/bin/env python3
"""Build and summarize a deterministic toy climate cube with Xarray."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def load_toolkit():
    try:
        import numpy as np
        import pandas as pd
        import xarray as xr
    except ImportError as exc:  # pragma: no cover - exercised via runtime tests
        raise SystemExit(
            "This script requires the scientific-python prefix at slurm/envs/scientific-python."
        ) from exc
    return np, pd, xr


def build_toy_cube():
    np, pd, xr = load_toolkit()
    time = pd.date_range("2025-01-01", periods=12, freq="MS")
    lat = np.array([-30.0, 0.0, 30.0], dtype=float)
    lon = np.array([0.0, 90.0, 180.0, 270.0], dtype=float)
    seasonal = 6.5 * np.sin(2 * np.pi * np.arange(12) / 12)
    lat_effect = np.array([2.0, 0.5, -1.5], dtype=float)
    lon_effect = np.array([0.0, 0.8, 1.1, 0.4], dtype=float)
    values = np.zeros((12, 3, 4), dtype=float)
    for i, seasonal_value in enumerate(seasonal):
        values[i, :, :] = 18.0 + seasonal_value + lat_effect[:, None] + lon_effect[None, :]
    return xr.DataArray(
        values,
        coords={"time": time, "lat": lat, "lon": lon},
        dims=("time", "lat", "lon"),
        name="temperature_c",
    )


def summarize_cube(cube) -> dict[str, Any]:
    np, _, _ = load_toolkit()
    monthly_mean = cube.mean(dim=("lat", "lon"))
    per_lat_mean = cube.mean(dim=("time", "lon"))
    hottest_index = int(np.argmax(monthly_mean.values))
    hottest_month = cube.coords["time"].values[hottest_index]
    hottest_month_str = str(hottest_month)[:7]
    return {
        "dims": {name: int(size) for name, size in cube.sizes.items()},
        "annual_mean_c": round(float(cube.mean().item()), 6),
        "seasonal_range_c": round(float(monthly_mean.max().item() - monthly_mean.min().item()), 6),
        "hottest_month": hottest_month_str,
        "per_lat_mean_c": [
            {"lat": round(float(lat), 3), "mean_c": round(float(value), 6)}
            for lat, value in zip(cube.coords["lat"].values, per_lat_mean.values, strict=True)
        ],
    }


def write_json(payload: dict[str, Any], out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()
    cube = build_toy_cube()
    payload = summarize_cube(cube)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
