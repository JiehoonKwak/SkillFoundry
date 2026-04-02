from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def load_metadata(skill_dir: Path) -> dict:
    return load_json(skill_dir / "metadata.yaml")


def default_example_path(skill_dir: Path) -> Path:
    return skill_dir / "examples" / "toy_input.json"


def write_tsv(path: Path, spec: dict) -> None:
    rows = spec.get("rows", [])
    columns = spec.get("columns")
    frame = pd.DataFrame(rows)
    if columns:
        for column in columns:
            if column not in frame.columns:
                frame[column] = ""
        frame = frame[columns]
    frame.to_csv(path, sep="\t", index=False)


def write_markdown(path: Path, spec: dict) -> None:
    title = spec.get("title") or path.stem.replace("_", " ").title()
    lines = [f"# {title}", ""]
    for section in spec.get("sections", []):
        name = section["name"]
        lines.append(f"## {name}")
        lines.append("")
        for item in section.get("bullets", []):
            lines.append(f"- {item}")
        for paragraph in section.get("paragraphs", []):
            lines.append(paragraph)
        lines.append("")
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def write_h5ad(path: Path, spec: dict) -> None:
    obs_names = spec.get("obs_names", ["obs_1", "obs_2", "obs_3"])
    var_names = spec.get("var_names", ["feature_1", "feature_2"])
    matrix = np.asarray(spec.get("matrix", [[1.0, 0.0], [0.0, 1.0], [0.5, 0.5]]), dtype=float)
    obs_frame = pd.DataFrame(spec.get("obs", {}), index=obs_names)
    var_frame = pd.DataFrame(spec.get("var", {}), index=var_names)
    adata = ad.AnnData(X=matrix, obs=obs_frame, var=var_frame)
    for key, values in spec.get("obsm", {}).items():
        adata.obsm[key] = np.asarray(values, dtype=float)
    adata.uns["pseudo_contract"] = {"generator": "experiments/sc_skills/shared/contract_runner.py"}
    adata.write_h5ad(path)


def write_h5mu(path: Path, spec: dict) -> None:
    with h5py.File(path, "w") as handle:
        handle.attrs["format"] = "pseudo_h5mu"
        handle.attrs["generator"] = "experiments/sc_skills/shared/contract_runner.py"
        modalities = handle.create_group("modalities")
        for modality_name, modality in spec.get("modalities", {}).items():
            group = modalities.create_group(modality_name)
            group.create_dataset("obs_names", data=np.asarray(modality.get("obs_names", []), dtype="S"))
            group.create_dataset("var_names", data=np.asarray(modality.get("var_names", []), dtype="S"))
            group.create_dataset("matrix", data=np.asarray(modality.get("matrix", []), dtype=float))
        obsm = handle.create_group("obsm")
        for key, values in spec.get("obsm", {}).items():
            obsm.create_dataset(key, data=np.asarray(values, dtype=float))


def write_deliverable(path: Path, spec: dict) -> None:
    kind = spec["kind"]
    path.parent.mkdir(parents=True, exist_ok=True)
    if kind == "tsv":
        write_tsv(path, spec)
        return
    if kind == "md":
        write_markdown(path, spec)
        return
    if kind == "h5ad":
        write_h5ad(path, spec)
        return
    if kind == "h5mu":
        write_h5mu(path, spec)
        return
    raise ValueError(f"Unsupported deliverable kind: {kind}")


def validate_outputs(skill_dir: Path, outdir: Path) -> None:
    metadata = load_metadata(skill_dir)
    for deliverable in metadata["deliverables"]:
        path = outdir / deliverable["path"]
        if not path.exists():
            raise AssertionError(f"Missing deliverable: {path}")
        if deliverable["kind"] == "tsv":
            frame = pd.read_csv(path, sep="\t")
            required_columns = deliverable.get("required_columns", [])
            missing = [column for column in required_columns if column not in frame.columns]
            if missing:
                raise AssertionError(f"Missing TSV columns in {path}: {missing}")
        elif deliverable["kind"] == "md":
            text = path.read_text(encoding="utf-8")
            for section in deliverable.get("required_sections", []):
                heading = f"## {section}"
                if heading not in text:
                    raise AssertionError(f"Missing markdown section {heading} in {path}")
        elif deliverable["kind"] == "h5ad":
            adata = ad.read_h5ad(path)
            if adata.n_obs == 0 or adata.n_vars == 0:
                raise AssertionError(f"Empty AnnData artifact: {path}")
        elif deliverable["kind"] == "h5mu":
            with h5py.File(path, "r") as handle:
                fmt = handle.attrs.get("format", "")
                if isinstance(fmt, bytes):
                    fmt = fmt.decode()
                if fmt != "pseudo_h5mu":
                    raise AssertionError(f"Unexpected pseudo h5mu marker in {path}")


def run_skill(skill_dir: Path, *, input_path: Path, outdir: Path) -> dict:
    payload = load_json(input_path)
    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    for deliverable_name, spec in payload["deliverables"].items():
        write_deliverable(outdir / deliverable_name, spec)
    validate_outputs(skill_dir, outdir)
    result = {
        "skill_dir": str(skill_dir),
        "input_path": str(input_path),
        "outdir": str(outdir),
        "written_files": sorted(str(path.relative_to(outdir)) for path in outdir.rglob("*") if path.is_file()),
    }
    (outdir / "run_summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main_for_skill(skill_dir: Path, argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=f"Run the toy contract for {skill_dir.name}.")
    parser.add_argument("--input", type=Path, default=default_example_path(skill_dir))
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    result = run_skill(skill_dir, input_path=args.input, outdir=args.outdir)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0
