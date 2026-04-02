#!/usr/bin/env python3
"""Generate a deterministic toy image and extract regionprops features."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def load_toolkit():
    try:
        import numpy as np
        from skimage.draw import ellipse
        from skimage.measure import label, regionprops_table
    except ImportError as exc:  # pragma: no cover - exercised via runtime tests
        raise SystemExit(
            "This script requires the scientific-python prefix at slurm/envs/scientific-python."
        ) from exc
    return np, ellipse, label, regionprops_table


def generate_toy_image(seed: int = 13, size: int = 112):
    np, ellipse, _, _ = load_toolkit()
    image = np.zeros((size, size), dtype=float)
    specs = [
        ((28, 30), 9, 15, 0.72),
        ((58, 74), 12, 10, 0.88),
        ((84, 40), 8, 13, 0.79),
    ]
    for center, r_radius, c_radius, intensity in specs:
        rr, cc = ellipse(center[0], center[1], r_radius, c_radius, shape=image.shape)
        image[rr, cc] = intensity
    rng = np.random.default_rng(seed)
    image += rng.normal(loc=0.0, scale=0.02, size=image.shape)
    return np.clip(image, 0.0, 1.0)


def extract_regionprops_features(image, threshold: float = 0.35) -> dict[str, Any]:
    np, _, label, regionprops_table = load_toolkit()
    mask = image > threshold
    labeled = label(mask)
    props = regionprops_table(
        labeled,
        intensity_image=image,
        properties=(
            "label",
            "area",
            "eccentricity",
            "perimeter",
            "axis_major_length",
            "axis_minor_length",
            "mean_intensity",
            "centroid",
        ),
    )
    objects = []
    for idx in range(len(props["label"])):
        objects.append(
            {
                "label": int(props["label"][idx]),
                "area": int(props["area"][idx]),
                "eccentricity": round(float(props["eccentricity"][idx]), 6),
                "perimeter": round(float(props["perimeter"][idx]), 6),
                "axis_major_length": round(float(props["axis_major_length"][idx]), 6),
                "axis_minor_length": round(float(props["axis_minor_length"][idx]), 6),
                "mean_intensity": round(float(props["mean_intensity"][idx]), 6),
                "centroid_row": round(float(props["centroid-0"][idx]), 6),
                "centroid_col": round(float(props["centroid-1"][idx]), 6),
            }
        )
    objects.sort(key=lambda item: (-item["area"], item["label"]))
    return {
        "image_shape": list(image.shape),
        "object_count": len(objects),
        "total_area": int(sum(item["area"] for item in objects)),
        "foreground_pixels": int(mask.sum()),
        "threshold": round(float(threshold), 6),
        "intensity_range": [round(float(np.min(image)), 6), round(float(np.max(image)), 6)],
        "mean_intensity_range": [
            round(float(min(item["mean_intensity"] for item in objects)), 6),
            round(float(max(item["mean_intensity"] for item in objects)), 6),
        ],
        "objects": objects,
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
    parser.add_argument("--seed", type=int, default=13, help="Random seed for deterministic noise.")
    parser.add_argument("--size", type=int, default=112, help="Toy image height/width in pixels.")
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.35,
        help="Foreground threshold applied before labeling.",
    )
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    if args.size <= 0:
        raise SystemExit("size must be a positive integer.")
    if not 0.0 < args.threshold < 1.0:
        raise SystemExit("threshold must be between 0 and 1.")

    image = generate_toy_image(seed=args.seed, size=args.size)
    payload = extract_regionprops_features(image, threshold=args.threshold)
    payload["seed"] = args.seed
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
