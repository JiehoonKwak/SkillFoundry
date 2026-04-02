#!/usr/bin/env python3
"""Generate a synthetic plant image and summarize a simple PlantCV segmentation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import imageio.v3 as iio
import numpy as np
from plantcv import plantcv as pcv


def build_synthetic_image() -> np.ndarray:
    image = np.zeros((64, 64, 3), dtype=np.uint8)
    image[16:48, 20:44, 1] = 180
    image[16:48, 20:44, 0] = 40
    image[16:48, 20:44, 2] = 20
    return image


def summarize_image(image: np.ndarray, threshold: int) -> tuple[dict[str, object], np.ndarray]:
    saturation = pcv.rgb2gray_hsv(rgb_img=image, channel="s")
    mask = pcv.threshold.binary(gray_img=saturation, threshold=threshold, object_type="light")
    ys, xs = np.where(mask > 0)
    if len(xs) == 0:
        raise SystemExit("Threshold removed the full foreground mask; choose a smaller threshold.")
    inside = saturation[mask > 0]
    payload = {
        "mode": "generated_toy",
        "image_shape": list(image.shape),
        "threshold": int(threshold),
        "foreground_pixel_count": int((mask > 0).sum()),
        "bbox": [int(xs.min()), int(ys.min()), int(xs.max()), int(ys.max())],
        "bbox_width": int(xs.max() - xs.min() + 1),
        "bbox_height": int(ys.max() - ys.min() + 1),
        "mean_saturation_inside": round(float(inside.mean()), 6),
        "max_mask_value": int(mask.max()),
    }
    return payload, mask.astype(np.uint8)


def write_json(payload: dict[str, object], out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def maybe_write_image(image: np.ndarray, out_path: Path | None) -> None:
    if out_path is None:
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    iio.imwrite(out_path, image)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threshold", type=int, default=10, help="Binary threshold in 0-255.")
    parser.add_argument("--image-out", type=Path, default=None, help="Optional PNG path for the synthetic RGB image.")
    parser.add_argument("--mask-out", type=Path, default=None, help="Optional PNG path for the segmentation mask.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON summary path.")
    args = parser.parse_args()

    if args.threshold < 0 or args.threshold > 255:
        raise SystemExit("--threshold must be between 0 and 255.")

    image = build_synthetic_image()
    payload, mask = summarize_image(image=image, threshold=args.threshold)
    maybe_write_image(image=image, out_path=args.image_out)
    maybe_write_image(image=mask, out_path=args.mask_out)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
