#!/usr/bin/env python3
"""Build a tiny deterministic EEG connectivity graph with MNE-Connectivity."""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]


def load_stack():
    try:
        import mne
        import numpy as np
        from mne_connectivity import spectral_connectivity_epochs
    except ImportError as exc:  # pragma: no cover - exercised through runtime invocation
        raise SystemExit(
            "This script requires slurm/envs/neuro. Run it with slurm/envs/neuro/bin/python."
        ) from exc
    return mne, np, spectral_connectivity_epochs


def build_epochs():
    mne, np, _ = load_stack()
    rng = np.random.default_rng(7)
    sfreq = 100.0
    n_epochs = 8
    n_times = 400
    t = np.arange(n_times) / sfreq
    base = np.sin(2 * np.pi * 10 * t)
    epochs_data = np.empty((n_epochs, 3, n_times), dtype=float)
    for idx in range(n_epochs):
        epochs_data[idx, 0, :] = base + 0.05 * rng.standard_normal(n_times)
        epochs_data[idx, 1, :] = 0.8 * base + 0.05 * rng.standard_normal(n_times)
        epochs_data[idx, 2, :] = np.sin(2 * np.pi * 20 * t + 0.3) + 0.05 * rng.standard_normal(n_times)
    info = mne.create_info(["Fz", "Cz", "Pz"], sfreq=sfreq, ch_types="eeg")
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="There were no Annotations stored.*", category=RuntimeWarning)
        epochs = mne.EpochsArray(epochs_data, info, verbose=False)
    return epochs


def summarize_graph(threshold: float) -> dict[str, object]:
    _, np, spectral_connectivity_epochs = load_stack()
    epochs = build_epochs()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="There were no Annotations stored.*", category=RuntimeWarning)
        conn = spectral_connectivity_epochs(
            epochs,
            method="coh",
            mode="multitaper",
            sfreq=float(epochs.info["sfreq"]),
            fmin=8,
            fmax=12,
            faverage=True,
            verbose=False,
        )
    matrix = conn.get_data(output="dense")[:, :, 0]
    channels = list(epochs.ch_names)
    rounded_matrix = [[round(float(value), 6) for value in row] for row in matrix]

    edges: list[dict[str, object]] = []
    max_edge: dict[str, object] | None = None
    for target_idx, target in enumerate(channels):
        for source_idx in range(target_idx):
            weight = round(float(matrix[target_idx, source_idx]), 6)
            source = channels[source_idx]
            if max_edge is None or weight > max_edge["weight"]:
                max_edge = {"source": source, "target": target, "weight": weight}
            if weight >= threshold:
                edges.append({"source": source, "target": target, "weight": weight})

    payload = {
        "channels": channels,
        "epochs": len(epochs),
        "samples_per_epoch": len(epochs.times),
        "frequency_band_hz": [8.0, 12.0],
        "threshold": threshold,
        "connectivity_matrix": rounded_matrix,
        "strong_edges": edges,
        "max_edge": max_edge,
    }
    return payload


def write_json(payload: dict[str, object], out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threshold", type=float, default=0.5, help="Minimum coherence to include in the edge list.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON summary path.")
    args = parser.parse_args()

    payload = summarize_graph(args.threshold)
    out_path = args.out if args.out is None or args.out.is_absolute() else ROOT / args.out
    write_json(payload, out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
