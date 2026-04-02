from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd


def _write_table(path: Path, spec: dict) -> None:
    frame = pd.DataFrame(spec["rows"], columns=spec["columns"])
    frame.to_csv(path, sep="\t", index=False)


def _write_markdown(path: Path, spec: dict) -> None:
    lines: list[str] = []
    title = spec.get("title")
    if title:
        lines.append(title)
        lines.append("")
    for section in spec.get("sections", []):
        lines.append(f"## {section['heading']}")
        lines.append("")
        lines.append(section["body"])
        lines.append("")
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def _write_h5ad(path: Path, spec: dict) -> None:
    obs = pd.DataFrame(index=spec["obs_names"])
    for column, values in spec.get("obs_columns", {}).items():
        obs[column] = values
    var = pd.DataFrame(index=spec["var_names"])
    for column, values in spec.get("var_columns", {}).items():
        var[column] = values
    adata = ad.AnnData(X=np.asarray(spec["X"], dtype=float), obs=obs, var=var)
    adata.write_h5ad(path)


def _write_h5mu(path: Path, spec: dict) -> None:
    with h5py.File(path, "w") as handle:
        handle.attrs["encoding-type"] = "MuData"
        handle.attrs["encoding-version"] = "0.1.0"
        handle.create_dataset("latent", data=np.asarray(spec["latent"], dtype=float))
        handle.create_dataset("obs_names", data=np.asarray(spec["obs_names"], dtype="S"))
        handle.create_dataset("latent_names", data=np.asarray(spec["latent_names"], dtype="S"))
        mods = handle.create_group("mod")
        for modality_name, modality in spec.get("modalities", {}).items():
            group = mods.create_group(modality_name)
            group.create_dataset("shape", data=np.asarray(modality["shape"], dtype=int))
            group.attrs["summary"] = modality["summary"]


def _validate_table(path: Path, spec: dict) -> None:
    frame = pd.read_csv(path, sep="\t")
    expected = spec["columns"]
    if list(frame.columns) != expected:
        raise AssertionError(f"{path} columns {list(frame.columns)} != {expected}")
    if len(frame) < spec.get("min_rows", 1):
        raise AssertionError(f"{path} has too few rows: {len(frame)}")


def _validate_markdown(path: Path, spec: dict) -> None:
    text = path.read_text(encoding="utf-8")
    title = spec.get("title")
    if title and title not in text:
        raise AssertionError(f"{path} missing title {title!r}")
    for section in spec.get("sections", []):
        heading = f"## {section['heading']}"
        if heading not in text:
            raise AssertionError(f"{path} missing heading {heading!r}")


def _validate_h5ad(path: Path, spec: dict) -> None:
    adata = ad.read_h5ad(path)
    expected_shape = (len(spec["obs_names"]), len(spec["var_names"]))
    if adata.shape != expected_shape:
        raise AssertionError(f"{path} shape {adata.shape} != {expected_shape}")


def _validate_h5mu(path: Path, spec: dict) -> None:
    with h5py.File(path, "r") as handle:
        if "latent" not in handle:
            raise AssertionError(f"{path} missing latent dataset")
        if "mod" not in handle:
            raise AssertionError(f"{path} missing modality group")
        latent_shape = tuple(handle["latent"].shape)
        expected = (len(spec["obs_names"]), len(spec["latent_names"]))
        if latent_shape != expected:
            raise AssertionError(f"{path} latent shape {latent_shape} != {expected}")


def load_example(task_root: Path) -> dict:
    return json.loads((task_root / "examples" / "toy_input.json").read_text(encoding="utf-8"))


def run_task(task_root: Path, output_dir: Path) -> None:
    spec = load_example(task_root)
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, table_spec in spec.get("tables", {}).items():
        _write_table(output_dir / name, table_spec)
    for name, markdown_spec in spec.get("markdown", {}).items():
        _write_markdown(output_dir / name, markdown_spec)
    for name, h5ad_spec in spec.get("h5ad", {}).items():
        _write_h5ad(output_dir / name, h5ad_spec)
    for name, h5mu_spec in spec.get("h5mu", {}).items():
        _write_h5mu(output_dir / name, h5mu_spec)


def validate_task(task_root: Path, output_dir: Path) -> None:
    spec = load_example(task_root)
    for name, table_spec in spec.get("tables", {}).items():
        _validate_table(output_dir / name, table_spec)
    for name, markdown_spec in spec.get("markdown", {}).items():
        _validate_markdown(output_dir / name, markdown_spec)
    for name, h5ad_spec in spec.get("h5ad", {}).items():
        _validate_h5ad(output_dir / name, h5ad_spec)
    for name, h5mu_spec in spec.get("h5mu", {}).items():
        _validate_h5mu(output_dir / name, h5mu_spec)


def cli(task_root: Path, mode: str) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if mode == "run":
        run_task(task_root, args.output_dir)
    else:
        validate_task(task_root, args.output_dir)
    return 0
