---
name: tissue-niche-annotation-portable-skill
description: Use this experiment skill to annotate tissue niches for MERFISH, Xenium, CosMx, and seqFISH style single-cell-resolution spatial data with a concrete runbook for dataset inspection, spatial graph construction, neighborhood profiling, candidate niche discovery, and missing-tool recovery.
---
## Purpose

Use this package when you already have single-cell-resolution spatial data with per-cell coordinates and at least one usable cell-type or state annotation column.

Do not use this package as the first step if cell-type labels are missing or unreliable. Return to the cell-type annotation workflow first. A common upstream input is `annotation_table.tsv` from the separate `annotation` skill, using `predicted_label` or `broad_label` as the niche-ready label column. Do not use this package for spot-based assays such as Visium, Slide-seq, or ST.

## What this package gives you

- A real-run execution playbook with concrete shell commands, Python snippets, intermediate artifacts, and decision rules.
- A deterministic offline surrogate runner in `scripts/run_tissue_niche_annotation.py` for contract validation and fallback rehearsal.
- A runtime checker in `scripts/check_runtime.py`.
- A run planner in `scripts/plan_niche_run.py` that scaffolds the expected files, graph parameters, and stop conditions before you touch the real data.

## Real-run directory contract

Create one run directory and keep every intermediate artifact there:

```bash
RUN_DIR="$PWD/scratch/tissue_niche_real_run"
mkdir -p "$RUN_DIR"/{qc,deliverables}
```

The playbook below expects these key files:

- `00_runtime_check.json`
- `01_dataset_profile.json`
- `02_niche_run_plan.json`
- `03_spatial_graph.tsv`
- `04_neighborhood_profiles.tsv`
- `05_niche_candidates.tsv`
- `06_niche_review.md`
- `niche_labels.tsv`
- `niche_markers.tsv`
- `tissue_niche_report.md`

## Step 0. Check the runtime before touching the data

Run:

```bash
python3 scripts/check_runtime.py --out "$RUN_DIR/qc/00_runtime_check.json"
```

Inspect:

- `conda_available`
- `missing_modules`
- `bootstrap_commands`
- `fallback_guidance`

Decision rules:

- If `squidpy` is missing, either install it or accept that graph construction and neighborhood enrichment must use the bundled surrogate path for a provisional run.
- If `cellcharter` is missing, you can still proceed with neighborhood-composition profiles and prototype-style niche scoring, but do not present the result as a CellCharter run.
- If `utag` is missing, use the graph-aware fallback path and mark niche boundaries as provisional.

Exact bootstrap suggestions:

Preferred if `conda` is available:

```bash
conda create -y -n sc-tissue-niche python=3.12
conda activate sc-tissue-niche
pip install --upgrade pip
pip install anndata scanpy pandas numpy scipy scikit-learn squidpy
```

Minimal venv route:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install anndata scanpy pandas numpy scipy scikit-learn squidpy
```

For CellCharter:

```bash
source .venv/bin/activate
pip install cellcharter
```

For UTAG:

```bash
python3 -m venv .venv-utag
source .venv-utag/bin/activate
pip install --upgrade pip
pip install 'scanpy>=1.9.3,<1.12' 'squidpy>=1.4'
pip install git+https://github.com/ElementoLab/utag.git@main
```

## Step 1. Inspect the annotated spatial dataset

Goal: confirm that the data are single-cell resolution, coordinates exist, and a niche-ready label column is available.

If the dataset came from the `annotation` skill, first check whether `predicted_label`, `broad_label`, `label_confidence`, and `marker_status` are present in `adata.obs` or can be merged from `annotation_table.tsv`.

Run:

```bash
python3 - <<'PY'
import json
from pathlib import Path

import anndata as ad

run_dir = Path("scratch/tissue_niche_real_run")
adata = ad.read_h5ad("query.h5ad", backed="r")
obs_cols = list(map(str, adata.obs.columns))
var_cols = list(map(str, adata.var.columns))
coord_pairs = [
    ("x_coord", "y_coord"),
    ("x_centroid", "y_centroid"),
    ("spatial_x", "spatial_y"),
    ("x", "y"),
]
coords = next(([x, y] for x, y in coord_pairs if x in obs_cols and y in obs_cols), [])
cell_type_candidates = [
    col for col in obs_cols
    if any(tok in col.lower() for tok in ["cell_type", "celltype", "broad_label", "predicted_label", "annotation", "state"])
]
sample_candidates = [
    col for col in obs_cols
    if any(tok in col.lower() for tok in ["sample", "batch", "patient", "donor", "fov", "roi", "slide"])
]
profile = {
    "adata_path": str(Path("query.h5ad").resolve()),
    "n_obs": int(adata.n_obs),
    "n_vars": int(adata.n_vars),
    "obs_columns": obs_cols,
    "var_columns": var_cols,
    "coordinate_keys": coords,
    "cell_type_candidates": cell_type_candidates,
    "sample_candidates": sample_candidates,
    "platform_hint": str(adata.uns.get("spatial_platform", adata.uns.get("platform", ""))),
}
run_dir.joinpath("qc/01_dataset_profile.json").write_text(json.dumps(profile, indent=2) + "\n", encoding="utf-8")
PY
```

Capture before moving on:

- `n_obs`, `n_vars`
- coordinate columns
- a label column such as `cell_type`, `broad_label`, or `predicted_label`
- a sample or slide column if multiple fields of view or slides are present
- platform hints in `uns` or file metadata

Stop conditions:

- stop if there is no coordinate pair
- stop if there is no usable cell-type/state label column
- stop if the platform is spot-based
- stop if the dataset mixes multiple slides but no slide/sample key is available
- stop if the only available labels are low-confidence or obviously provisional; return to the annotation workflow first

## Step 2. Scaffold the niche run plan before graph construction

Use the helper script so the graph parameters and expected artifacts are explicit:

```bash
python3 scripts/plan_niche_run.py \
  --adata query.h5ad \
  --cell-type-column cell_type \
  --x-key x_coord \
  --y-key y_coord \
  --sample-key sample_id \
  --graph-mode knn \
  --k 12 \
  --out "$RUN_DIR/qc/02_niche_run_plan.json"
```

Inspect:

- chosen graph mode
- `k` or `radius`
- expected outputs
- declared stop conditions

Decision rules:

- default to `kNN` if density varies between regions or if coordinates are in arbitrary units
- use `radius` only when coordinates are already calibrated and a biologically meaningful distance threshold is known
- for very small panels or tiny regions, reduce `k` so each cell does not connect to nearly the entire dataset

## Step 3. Build the spatial graph

Default public path with Squidpy:

```python
import scanpy as sc
import squidpy as sq

adata = sc.read_h5ad("query.h5ad")
sq.gr.spatial_neighbors(adata, coord_type="generic", spatial_key="spatial", n_neighs=12)
adata.write_h5ad("scratch/tissue_niche_real_run/qc/03_graph_ready.h5ad")
```

If the dataset stores coordinates in `obs`, first copy them:

```python
import numpy as np

adata.obsm["spatial"] = adata.obs[["x_coord", "y_coord"]].to_numpy(dtype=float)
```

Expected artifact:

- a graph-ready H5AD or a TSV export such as `03_spatial_graph.tsv`

If `squidpy` is missing:

- use the bundled deterministic graph logic from `scripts/run_tissue_niche_annotation.py` as a fallback rehearsal path
- document that the fallback validates workflow shape but does not replace Squidpy neighborhood enrichment on real data

## Step 4. Compute neighborhood composition and enrichment

Default public path:

```python
import pandas as pd
import squidpy as sq

sq.gr.nhood_enrichment(adata, cluster_key="cell_type")
sq.gr.co_occurrence(adata, cluster_key="cell_type")
```

Also build explicit per-cell neighborhood composition features, because they are easier to interpret and easier to use in fallback mode. Save:

- `04_neighborhood_profiles.tsv`
  - per-cell fractions or weighted counts of nearby cell types
- `04_nhood_enrichment.tsv`
  - cell-type pair enrichment summaries when Squidpy is available

Decision rules:

- use raw local composition when niche boundaries look sharp
- use graph-smoothed composition when local neighborhoods are noisy or sparsely sampled
- keep boundary scores or entropy-like metrics explicit instead of collapsing them into a hard niche label too early

## Step 5. Discover candidate niches

Preferred public paths:

- CellCharter: neighborhood aggregation plus clustering
- UTAG: graph-aware microenvironment grouping

Representative CellCharter-style path:

```python
# pseudocode-level public path
# 1. aggregate neighborhood features
# 2. cluster the aggregated features
# 3. review cluster size and composition
```

Representative UTAG-style path:

```python
# prepare AnnData with obsm['spatial']
# run UTAG according to the pinned repo instructions
# review cluster labels and spatial coherence
```

If `cellcharter` and `utag` are both unavailable:

- use the bundled fallback path already implemented in `scripts/run_tissue_niche_annotation.py`
- keep the result labeled as a surrogate or provisional niche call
- preserve intermediate outputs such as smoothed neighborhood profiles and prototype scores

Expected artifact:

- `05_niche_candidates.tsv`
  - per-cell candidate niche labels
  - confidence or margin
  - boundary/interface score
  - dominant neighbor type

## Step 6. Interpret and annotate the niches

For each candidate niche, answer:

- which cell types dominate locally?
- is the niche coherent within one region or spread across multiple slides?
- is it an interface niche or a core niche?
- do marker or state summaries support the interpretation?

Write:

- `06_niche_review.md`
- `niche_markers.tsv`
- `tissue_niche_report.md`

Do not simply rename `cell_type` values as niches. A niche must reflect neighborhood context, not just the source label of the focal cell.

## Step 7. Final deliverables

Required final outputs:

- `niche_labels.tsv`
- `niche_markers.tsv`
- `tissue_niche_report.md`

Minimum required columns in `niche_labels.tsv`:

- `cell_id`
- `sample_id`
- `cell_type`
- `niche_label`
- `niche_confidence`
- `dominant_neighbor_type`
- `boundary_score`
- `notes`

## Missing-tool fallback rules

- If `squidpy` is missing:
  - install it if possible
  - otherwise use the bundled surrogate graph path for provisional niche exploration
- If `cellcharter` is missing:
  - keep the run on neighborhood profiles + explicit prototype or heuristic scoring
  - mark the report as “no CellCharter clustering available”
- If `utag` is missing:
  - keep the graph-aware fallback path
  - record that the run used surrogate niche grouping rather than UTAG proper
- If the label column is weak:
  - return to the annotation workflow and do not continue with niche labeling
- Do not repeat CELLxGENE reference search or label transfer in this package. Those belong in the separate `annotation` skill.

## Toy validation

The toy path remains for deterministic validation only:

```bash
python3 scripts/run_tissue_niche_annotation.py --input examples/toy_input.json --outdir scratch/tissue-niche-annotation_toy_run
python3 scripts/validate_outputs.py --outdir scratch/tissue-niche-annotation_toy_run
python3 -m unittest tests/test_contract.py tests/test_tissue_niche_annotation.py
```

Treat these outputs as contract validation, not biological evidence.

## Starter scope

- Truly computed: spatial kNN graph construction, neighborhood composition, one-step graph smoothing, prototype-based niche scoring, and niche-vs-rest marker summaries on small synthetic inputs.
- Approximated: Squidpy neighborhood enrichment, CellCharter aggregation/clustering, and UTAG grouping are represented by documented public paths plus deterministic bundled fallbacks.
- Surrogate only: the toy niche prototypes and synthetic marker values validate workflow shape; they do not replace dataset-specific calibration or biological review.

## Hand-off contract

Before finishing, leave behind:

- the declared deliverables
- the QC files listed above
- the run planner JSON
- enough provenance in `tissue_niche_report.md` for another shell-capable agent to rerun the workflow
