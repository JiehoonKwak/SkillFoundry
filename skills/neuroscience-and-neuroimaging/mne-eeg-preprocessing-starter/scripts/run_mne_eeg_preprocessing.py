#!/usr/bin/env python3
"""Create a synthetic EEG RawArray and summarize a simple MNE filter pass."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mne
import numpy as np


def build_raw() -> mne.io.RawArray:
    sfreq = 100.0
    time_axis = np.arange(400) / sfreq
    data = np.vstack(
        [
            np.sin(2 * np.pi * 10 * time_axis) + 0.4 * np.sin(2 * np.pi * 0.2 * time_axis),
            0.8 * np.cos(2 * np.pi * 12 * time_axis) + 0.3 * np.sin(2 * np.pi * 0.2 * time_axis + 0.3),
        ]
    ) * 1e-6
    info = mne.create_info(["Fz", "Cz"], sfreq=sfreq, ch_types=["eeg", "eeg"])
    return mne.io.RawArray(data, info, verbose="error")


def build_summary() -> dict[str, object]:
    raw = build_raw()
    filtered = raw.copy().filter(l_freq=1.0, h_freq=30.0, fir_design="firwin", verbose="error")
    before = raw.get_data()
    after = filtered.get_data()
    ratios = [round(float(a.std() / b.std()), 6) for a, b in zip(after, before)]
    return {
        "channel_names": raw.ch_names,
        "channel_count": len(raw.ch_names),
        "sfreq": round(float(raw.info["sfreq"]), 6),
        "duration_s": round(float(raw.times[-1] + 1.0 / raw.info["sfreq"]), 6),
        "pre_means": [round(float(value), 12) for value in before.mean(axis=1)],
        "post_means": [round(float(value), 12) for value in after.mean(axis=1)],
        "filtered_std_ratio": ratios,
        "first_channel_first_samples": [round(float(value), 12) for value in after[0, :5]],
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

    payload = build_summary()
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
