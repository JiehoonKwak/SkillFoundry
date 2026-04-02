#!/usr/bin/env python3
"""Summarize a deterministic toy posterior with ArviZ diagnostics."""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path


def load_toolkit():
    warnings.filterwarnings(
        "ignore",
        message=".*ArviZ is undergoing a major refactor.*",
        category=FutureWarning,
    )
    try:
        import arviz as az
        import numpy as np
    except ImportError as exc:  # pragma: no cover - exercised via runtime tests
        raise SystemExit(
            "This script requires the dedicated statistics prefix at slurm/envs/statistics."
        ) from exc
    return az, np


def build_toy_inferencedata(draws: int, chains: int, seed: int):
    az, np = load_toolkit()
    rng = np.random.default_rng(seed)
    alpha = rng.normal(loc=1.2, scale=0.05, size=(chains, draws))
    beta = rng.normal(loc=-0.75, scale=0.08, size=(chains, draws))
    sigma = np.abs(rng.normal(loc=0.35, scale=0.025, size=(chains, draws)))
    return az.from_dict(posterior={"alpha": alpha, "beta": beta, "sigma": sigma})


def summarize_diagnostics(idata) -> dict:
    az, _ = load_toolkit()
    summary = az.summary(idata, var_names=["alpha", "beta", "sigma"], hdi_prob=0.9, round_to=6)
    parameters = {}
    for name in ["alpha", "beta", "sigma"]:
        parameters[name] = {
            "mean": round(float(summary.loc[name, "mean"]), 6),
            "sd": round(float(summary.loc[name, "sd"]), 6),
            "hdi_90": [
                round(float(summary.loc[name, "hdi_5%"]), 6),
                round(float(summary.loc[name, "hdi_95%"]), 6),
            ],
            "ess_bulk": round(float(summary.loc[name, "ess_bulk"]), 6),
            "ess_tail": round(float(summary.loc[name, "ess_tail"]), 6),
            "r_hat": round(float(summary.loc[name, "r_hat"]), 6),
        }
    return {
        "parameter_count": len(parameters),
        "parameters": parameters,
        "max_rhat": round(max(item["r_hat"] for item in parameters.values()), 6),
        "min_ess_bulk": round(min(item["ess_bulk"] for item in parameters.values()), 6),
        "min_ess_tail": round(min(item["ess_tail"] for item in parameters.values()), 6),
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
    parser.add_argument("--draws", type=int, default=300, help="Posterior draws per chain.")
    parser.add_argument("--chains", type=int, default=4, help="Number of chains.")
    parser.add_argument("--seed", type=int, default=20260312, help="Deterministic random seed.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    if args.draws <= 0 or args.chains < 2:
        raise SystemExit("draws must be positive and chains must be at least 2.")

    idata = build_toy_inferencedata(draws=args.draws, chains=args.chains, seed=args.seed)
    payload = summarize_diagnostics(idata)
    payload["draws"] = args.draws
    payload["chains"] = args.chains
    payload["seed"] = args.seed
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
