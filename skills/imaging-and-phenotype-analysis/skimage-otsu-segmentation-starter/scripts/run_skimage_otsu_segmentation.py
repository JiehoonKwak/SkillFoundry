#!/usr/bin/env python3
"""Generate a toy image and summarize an Otsu-threshold segmentation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def load_toolkit():
    try:
        import numpy as np
        from skimage.draw import disk
        from skimage.filters import gaussian, threshold_otsu
        from skimage.measure import label, regionprops_table
    except ImportError as exc:  # pragma: no cover - exercised via runtime tests
        raise SystemExit(
            "This script requires the scientific-python prefix at slurm/envs/scientific-python."
        ) from exc
    return np, disk, gaussian, threshold_otsu, label, regionprops_table


def generate_toy_image(seed: int = 7, size: int = 96):
    np, disk, *_ = load_toolkit()
    image = np.zeros((size, size), dtype=float)
    disks = [
        ((24, 28), 10, 0.82),
        ((54, 62), 12, 0.92),
        ((72, 26), 9, 0.78),
    ]
    for center, radius, intensity in disks:
        rr, cc = disk(center, radius, shape=image.shape)
        image[rr, cc] = intensity
    rng = np.random.default_rng(seed)
    image += rng.normal(loc=0.0, scale=0.035, size=image.shape)
    return np.clip(image, 0.0, 1.0)


def segment_image(image) -> dict[str, Any]:
    np, _, gaussian, threshold_otsu, label, regionprops_table = load_toolkit()
    blurred = gaussian(image, sigma=1.0, preserve_range=True)
    threshold = float(threshold_otsu(blurred))
    initial_labels = label(blurred > threshold)
    initial_props = regionprops_table(initial_labels, properties=("label", "area"))
    keep_labels = [
        int(initial_props["label"][idx])
        for idx in range(len(initial_props["label"]))
        if int(initial_props["area"][idx]) >= 80
    ]
    mask = np.isin(initial_labels, keep_labels)
    labeled = label(mask)
    props = regionprops_table(labeled, properties=("label", "area", "centroid"))
    objects = []
    for idx in range(len(props["label"])):
        objects.append(
            {
                "label": int(props["label"][idx]),
                "area": int(props["area"][idx]),
                "centroid_row": round(float(props["centroid-0"][idx]), 4),
                "centroid_col": round(float(props["centroid-1"][idx]), 4),
            }
        )
    return {
        "image_shape": list(image.shape),
        "intensity_range": [round(float(np.min(image)), 6), round(float(np.max(image)), 6)],
        "threshold": round(threshold, 6),
        "object_count": len(objects),
        "foreground_pixels": int(mask.sum()),
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
    parser.add_argument("--seed", type=int, default=7, help="Random seed for toy image generation.")
    parser.add_argument("--size", type=int, default=96, help="Toy image height/width in pixels.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    if args.size <= 0:
        raise SystemExit("size must be a positive integer.")

    image = generate_toy_image(seed=args.seed, size=args.size)
    payload = segment_image(image)
    payload["seed"] = args.seed
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
