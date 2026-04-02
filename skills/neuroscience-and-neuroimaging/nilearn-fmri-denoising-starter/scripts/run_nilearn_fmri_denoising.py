#!/usr/bin/env python3
"""Create toy fMRI-like signals and summarize Nilearn confound cleaning."""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
from nilearn.signal import clean


def build_toy_timeseries() -> tuple[np.ndarray, np.ndarray]:
    frames = 60
    time_axis = np.linspace(0.0, 1.0, frames)
    confounds = np.column_stack(
        [
            np.linspace(-1.0, 1.0, frames),
            np.sin(2 * np.pi * time_axis),
        ]
    )
    signals = np.column_stack(
        [
            0.7 * confounds[:, 0] + np.sin(2 * np.pi * 3 * time_axis),
            -0.5 * confounds[:, 1] + np.cos(2 * np.pi * 4 * time_axis),
        ]
    )
    return signals, confounds


def summarize_cleaning() -> dict[str, object]:
    signals, confounds = build_toy_timeseries()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        cleaned = clean(
            signals,
            confounds=confounds,
            detrend=True,
            standardize="zscore_sample",
        )

    before_values = []
    after_values = []
    for column in range(signals.shape[1]):
        before_values.extend(abs(np.corrcoef(signals[:, column], confounds[:, idx])[0, 1]) for idx in range(confounds.shape[1]))
        after_values.extend(abs(np.corrcoef(cleaned[:, column], confounds[:, idx])[0, 1]) for idx in range(confounds.shape[1]))

    return {
        "frame_count": int(signals.shape[0]),
        "voxel_count": int(signals.shape[1]),
        "pre_confound_abs_correlation_mean": round(float(np.mean(before_values)), 6),
        "post_confound_abs_correlation_mean": round(float(np.mean(after_values)), 6),
        "post_confound_abs_correlation_max": round(float(np.max(after_values)), 6),
        "cleaned_means": [round(float(value), 6) for value in cleaned.mean(axis=0)],
        "cleaned_std_sample": [round(float(value), 6) for value in cleaned.std(axis=0, ddof=1)],
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
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    payload = summarize_cleaning()
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
