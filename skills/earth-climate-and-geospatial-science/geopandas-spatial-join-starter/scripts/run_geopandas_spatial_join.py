#!/usr/bin/env python3
"""Run a deterministic GeoPandas spatial join on a toy point-and-region dataset."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path


def configure_geospatial_env(prefix: Path | None = None) -> Path:
    if prefix is None:
        prefix = Path(sys.executable).resolve().parents[1]
    proj_dir = prefix / "share" / "proj"
    gdal_dir = prefix / "share" / "gdal"
    if proj_dir.exists():
        os.environ["PROJ_LIB"] = str(proj_dir)
        os.environ["PROJ_DATA"] = str(proj_dir)
    if gdal_dir.exists():
        os.environ["GDAL_DATA"] = str(gdal_dir)
    return prefix


def load_geospatial_stack():
    configure_geospatial_env()
    try:
        import geopandas as gpd
        from shapely.geometry import Point, box
    except ImportError as exc:  # pragma: no cover - exercised via runtime tests
        raise SystemExit(
            "This script requires the dedicated geospatial prefix at slurm/envs/geospatial. "
            "Run it with slurm/envs/geospatial/bin/python."
        ) from exc
    return gpd, Point, box


def build_regions():
    gpd, _, box = load_geospatial_stack()
    records = [
        {"region_id": "north", "region_name": "North Basin", "geometry": box(-1.0, 0.0, 2.0, 3.0)},
        {"region_id": "south", "region_name": "South Basin", "geometry": box(-1.0, -3.0, 2.0, 0.0)},
    ]
    return gpd.GeoDataFrame(records, crs="EPSG:4326")


def build_points():
    gpd, Point, _ = load_geospatial_stack()
    records = [
        {"point_id": "P1", "measurement": 10.5, "geometry": Point(0.4, 1.6)},
        {"point_id": "P2", "measurement": 9.2, "geometry": Point(1.4, 0.8)},
        {"point_id": "P3", "measurement": 7.8, "geometry": Point(0.1, -1.3)},
        {"point_id": "P4", "measurement": 5.0, "geometry": Point(3.5, 0.2)},
    ]
    return gpd.GeoDataFrame(records, crs="EPSG:4326")


def summarize_join(points, regions) -> dict:
    gpd, _, _ = load_geospatial_stack()
    joined = gpd.sjoin(
        points,
        regions[["region_id", "region_name", "geometry"]],
        how="left",
        predicate="within",
    ).sort_values("point_id")
    projected = joined.to_crs(epsg=3857)
    assigned_counts = {
        region_id: int((joined["region_id"] == region_id).sum())
        for region_id in regions["region_id"].tolist()
    }
    measurement_by_region = {
        region_id: round(
            float(joined.loc[joined["region_id"] == region_id, "measurement"].sum()),
            4,
        )
        for region_id in regions["region_id"].tolist()
    }
    joined_rows = []
    for row in joined.itertuples(index=False):
        joined_rows.append(
            {
                "point_id": row.point_id,
                "measurement": round(float(row.measurement), 4),
                "region_id": row.region_id,
                "region_name": row.region_name,
            }
        )
    return {
        "point_count": int(points.shape[0]),
        "region_count": int(regions.shape[0]),
        "assigned_counts": assigned_counts,
        "measurement_by_region": measurement_by_region,
        "unassigned_point_ids": joined.loc[joined["region_id"].isna(), "point_id"].tolist(),
        "projected_crs": str(projected.crs),
        "projected_total_bounds": [round(float(value), 3) for value in projected.total_bounds],
        "environment": {
            "proj_lib": os.environ.get("PROJ_LIB"),
            "gdal_data": os.environ.get("GDAL_DATA"),
        },
        "joined_rows": joined_rows,
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
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    regions = build_regions()
    points = build_points()
    payload = summarize_join(points, regions)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
