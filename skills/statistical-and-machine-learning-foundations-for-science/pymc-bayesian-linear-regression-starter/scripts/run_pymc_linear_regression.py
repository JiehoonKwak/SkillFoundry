#!/usr/bin/env python3
"""Fit a small Bayesian linear regression with PyMC on deterministic toy data."""

from __future__ import annotations

import argparse
import csv
import json
import warnings
from pathlib import Path


DEFAULT_INPUT = Path(__file__).resolve().parents[1] / "examples" / "toy_observations.tsv"


def load_modeling_stack():
    warnings.filterwarnings(
        "ignore",
        message=".*ArviZ is undergoing a major refactor.*",
        category=FutureWarning,
    )
    try:
        import arviz as az
        import numpy as np
        import pymc as pm
    except ImportError as exc:  # pragma: no cover - exercised via runtime tests
        raise SystemExit(
            "This script requires the dedicated statistics prefix at slurm/envs/statistics. "
            "Run it with slurm/envs/statistics/bin/python."
        ) from exc
    return az, np, pm


def load_observations(path: Path) -> tuple[list[float], list[float]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"No observations found in {path}")
    x_values = [float(row["x"]) for row in rows]
    y_values = [float(row["y"]) for row in rows]
    return x_values, y_values


def fit_linear_regression(
    x_values: list[float],
    y_values: list[float],
    draws: int,
    tune: int,
    chains: int,
    seed: int,
):
    az, np, pm = load_modeling_stack()
    x_array = np.asarray(x_values, dtype=float)
    y_array = np.asarray(y_values, dtype=float)
    with pm.Model() as model:
        intercept = pm.Normal("intercept", mu=0.0, sigma=5.0)
        slope = pm.Normal("slope", mu=0.0, sigma=5.0)
        sigma = pm.HalfNormal("sigma", sigma=1.0)
        mu = intercept + slope * x_array
        pm.Normal("obs", mu=mu, sigma=sigma, observed=y_array)
        idata = pm.sample(
            draws=draws,
            tune=tune,
            chains=chains,
            cores=1,
            random_seed=seed,
            init="adapt_diag",
            progressbar=False,
            target_accept=0.9,
            compute_convergence_checks=False,
            discard_tuned_samples=True,
        )
    summary = az.summary(idata, var_names=["intercept", "slope", "sigma"], hdi_prob=0.9, round_to=6)
    return idata, summary


def build_summary(idata, summary_frame, x_values: list[float], y_values: list[float], draws: int, tune: int, chains: int) -> dict:
    _, np, _ = load_modeling_stack()
    posterior_means = {
        name: round(float(summary_frame.loc[name, "mean"]), 6)
        for name in ["intercept", "slope", "sigma"]
    }
    hdi_90 = {
        name: [
            round(float(summary_frame.loc[name, "hdi_5%"]), 6),
            round(float(summary_frame.loc[name, "hdi_95%"]), 6),
        ]
        for name in ["intercept", "slope", "sigma"]
    }
    posterior_x = np.asarray([0.0, 2.5, 5.0], dtype=float)
    posterior_y = posterior_means["intercept"] + posterior_means["slope"] * posterior_x
    return {
        "observation_count": len(x_values),
        "draws": draws,
        "tune": tune,
        "chains": chains,
        "posterior_means": posterior_means,
        "hdi_90": hdi_90,
        "posterior_predictive_mean": [
            {"x": round(float(x_value), 3), "mean_y": round(float(y_value), 6)}
            for x_value, y_value in zip(posterior_x, posterior_y, strict=True)
        ],
        "observed_range": {
            "x": [min(x_values), max(x_values)],
            "y": [min(y_values), max(y_values)],
        },
        "posterior_sample_count": int(idata.posterior.sizes["chain"] * idata.posterior.sizes["draw"]),
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
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Tab-separated x/y observations.")
    parser.add_argument("--draws", type=int, default=150, help="Posterior draws per chain.")
    parser.add_argument("--tune", type=int, default=150, help="Tuning steps per chain.")
    parser.add_argument("--chains", type=int, default=2, help="Number of chains.")
    parser.add_argument("--seed", type=int, default=20260312, help="Random seed for deterministic sampling.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    if args.draws <= 0 or args.tune <= 0 or args.chains <= 0:
        raise SystemExit("draws, tune, and chains must all be positive integers.")

    x_values, y_values = load_observations(args.input)
    idata, summary_frame = fit_linear_regression(
        x_values=x_values,
        y_values=y_values,
        draws=args.draws,
        tune=args.tune,
        chains=args.chains,
        seed=args.seed,
    )
    payload = build_summary(
        idata=idata,
        summary_frame=summary_frame,
        x_values=x_values,
        y_values=y_values,
        draws=args.draws,
        tune=args.tune,
        chains=args.chains,
    )
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
