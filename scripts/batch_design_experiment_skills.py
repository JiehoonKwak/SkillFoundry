#!/usr/bin/env python3
"""Materialize and harden portable experiment skills under experiments/sc_skills/."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPERIMENT_ROOT = ROOT / "experiments" / "sc_skills"
MANIFEST_PATH = EXPERIMENT_ROOT / "batch_design_manifest.json"


RESOURCE_GROUP_ADDITIONS = {
    "bioinformatics_utilities": [
        {
            "label": "Biopython tutorial",
            "url": "https://biopython.org/docs/latest/Tutorial/index.html",
            "note": "Portable sequence parsing, feature extraction, and record-normalization surface for general coding agents.",
        },
        {
            "label": "HGNC REST API",
            "url": "https://www.genenames.org/help/rest/",
            "note": "Gene symbol normalization surface for human gene-centric lookup and panel curation.",
        },
    ],
    "communication": [
        {
            "label": "LIANA+ paper",
            "url": "https://www.nature.com/articles/s41556-024-01469-w",
            "note": "Primary paper describing the LIANA+ unified framework for cell-cell communication inference.",
        },
        {
            "label": "OmniPath paper",
            "url": "https://www.embopress.org/doi/full/10.15252/msb.20209923",
            "note": "Primary paper for the integrated intra- and intercellular signaling knowledge base used by LIANA and Squidpy communication workflows.",
        },
    ],
    "spatial_foundation": [
        {
            "label": "SpatialData documentation",
            "url": "https://spatialdata.scverse.org/en/stable/",
            "note": "scverse data model and IO layer for portable spatial single-cell analysis workflows.",
        },
        {
            "label": "Tangram paper",
            "url": "https://www.nature.com/articles/s41592-021-01264-7",
            "note": "Primary paper for single-cell to spatial alignment, projection, and deconvolution with Tangram.",
        },
        {
            "label": "CellTypist documentation",
            "url": "https://celltypist.readthedocs.io/en/stable/",
            "note": "Official usage guide for portable query-to-reference cell-type annotation runs.",
        },
    ],
    "spatial_structure": [
        {
            "label": "stLearn documentation",
            "url": "https://stlearn.readthedocs.io/en/latest/",
            "note": "Additional spatial statistics and neighborhood-analysis surface for domain and microenvironment workflows.",
        },
        {
            "label": "CellCharter documentation",
            "url": "https://cellcharter.readthedocs.io/en/latest/",
            "note": "Reference documentation for spatial clustering and cellular niche identification from single-cell-resolution spatial data.",
        },
        {
            "label": "Single-cell best practices neighborhood analysis",
            "url": "https://www.sc-best-practices.org/spatial/neighborhood.html",
            "note": "Public best-practices tutorial covering spatial graphs, neighborhood enrichment, and niche-oriented downstream interpretation.",
        },
        {
            "label": "Squidpy paper",
            "url": "https://www.nature.com/articles/s41592-021-01358-2",
            "note": "Primary paper for scalable spatial graph, neighborhood, and co-occurrence analysis.",
        },
        {
            "label": "CellCharter paper",
            "url": "https://doi.org/10.1038/s41587-023-01844-y",
            "note": "Primary paper for spatial clustering and characterization of tissue niches across spatial -omics data.",
        },
        {
            "label": "SpaGCN paper",
            "url": "https://www.nature.com/articles/s41592-021-01255-8",
            "note": "Primary graph-convolution method paper for spatial domain detection and domain-guided markers.",
        },
        {
            "label": "GraphST paper",
            "url": "https://www.nature.com/articles/s41467-023-36796-3",
            "note": "Primary paper for spatially informed clustering, integration, and scRNA-to-spatial transfer with GraphST.",
        },
    ],
    "multimodal_modeling": [
        {
            "label": "scArches documentation",
            "url": "https://scarches.readthedocs.io/en/latest/",
            "note": "Transfer-learning and reference-mapping extension to scvi-tools for multimodal modeling workflows.",
        },
        {
            "label": "MultiVI paper",
            "url": "https://www.nature.com/articles/s41592-023-01909-9",
            "note": "Primary paper for multimodal single-cell integration and uncertainty-aware latent modeling in scvi-tools.",
        },
        {
            "label": "scArches paper",
            "url": "https://www.nature.com/articles/s41587-021-01001-7",
            "note": "Primary paper for transfer-learning based mapping to single-cell reference atlases.",
        },
        {
            "label": "Single-cell multi-omics benchmark",
            "url": "https://www.nature.com/articles/s41592-024-02429-w",
            "note": "Benchmark comparing multimodal integration and prediction methods, useful for method selection and validation planning.",
        },
    ],
    "spatial_deconvolution": [
        {
            "label": "cell2location paper",
            "url": "https://www.nature.com/articles/s41587-021-01139-4",
            "note": "Primary paper for probabilistic cell-type abundance estimation in spatial transcriptomics.",
        },
        {
            "label": "DestVI paper",
            "url": "https://www.nature.com/articles/s41587-022-01272-8",
            "note": "Primary paper for modeling continuous cell states in spatial transcriptomics with scvi-tools.",
        },
        {
            "label": "Spatial integration benchmark",
            "url": "https://www.nature.com/articles/s41592-022-01480-9",
            "note": "Benchmark paper comparing mapping, transcript distribution prediction, and deconvolution methods.",
        },
    ],
    "tangram_extensions": [
        {
            "label": "Tangram paper",
            "url": "https://www.nature.com/articles/s41592-021-01264-7",
            "note": "Primary paper for Tangram-based mapping, projection, and dropout-aware expression transfer.",
        },
        {
            "label": "Spatial integration benchmark",
            "url": "https://www.nature.com/articles/s41592-022-01480-9",
            "note": "Benchmark paper for selecting and sanity-checking Tangram-style extensions against alternative methods.",
        },
    ],
    "trajectory_modeling": [
        {
            "label": "scVelo paper",
            "url": "https://www.nature.com/articles/s41587-020-0591-3",
            "note": "Primary paper for dynamical RNA velocity modeling and transient-state trajectory inference.",
        },
        {
            "label": "CellRank paper",
            "url": "https://www.nature.com/articles/s41592-021-01346-6",
            "note": "Primary paper for directed cell-fate mapping and terminal-state probability estimation.",
        },
    ],
}


FAMILY_WORKFLOW = {
    "spatial_foundation": [
        "Inspect the spatial object, coordinate fields, metadata columns, and gene overlap with any reference before attempting transfer.",
        "Choose an explicit public reference source and keep the choice reproducible.",
        "Perform label transfer or marker-guided annotation with deterministic saved inputs and outputs.",
        "Validate labels with markers and spatial neighborhood context before final interpretation.",
    ],
    "tangram_extensions": [
        "Confirm reference-to-spatial feature compatibility before mapping or projection.",
        "Use Tangram-style mapping or projection with explicit QC and saved parameters.",
        "Measure held-out agreement, marker recovery, or assignment confidence on toy outputs.",
        "Document exactly how the extension differs from plain mapping.",
    ],
    "spatial_structure": [
        "Build the spatial graph from explicit coordinates and sample boundaries.",
        "Run domain or neighborhood statistics with transparent parameters.",
        "Write machine-readable labels and markers alongside the narrative report.",
        "Treat low-support domains as exploratory rather than final biology.",
    ],
    "spatial_deconvolution": [
        "Profile the spot-level object and the single-cell reference before choosing a model.",
        "Keep the abundance or assignment matrix as a first-class artifact.",
        "Compare at least one alternate method or sanity-check summary.",
        "Document how pseudo-data differs from a real abundance inference run.",
    ],
    "communication": [
        "Normalize group labels and metadata before running any ligand-receptor analysis.",
        "Keep ranked interaction tables and prioritized summaries separate.",
        "Record the evidence source or scoring method for every prioritized pair.",
        "Flag spatial support as optional rather than guaranteed.",
    ],
    "multimodal_modeling": [
        "Ensure observation identifiers align across modalities before integration.",
        "Materialize a latent representation plus modality-specific QC outputs.",
        "Keep modality provenance explicit in the final report.",
        "Treat toy multimodal outputs as contract validation, not biological discovery.",
    ],
    "trajectory_modeling": [
        "Prepare coordinates, neighborhood structure, or latent inputs before pseudotime or fate estimation.",
        "Write machine-readable trajectory coordinates and fate summaries.",
        "Keep uncertainty or caveat notes explicit in the report.",
        "Use toy lineage progressions only for contract validation.",
    ],
    "bioinformatics_utilities": [
        "Normalize identifiers before cross-source lookup or reconciliation.",
        "Prefer public APIs and documented CLIs over hidden wrappers.",
        "Return structured tables and provenance notes, not only prose.",
        "Document ambiguity, fallback behavior, and unresolved aliases.",
    ],
}


COMMON_PROMPT_RULES = [
    "Do not use `examples/toy_input.json` to store final deliverables or narrative report text. The toy input may only hold raw synthetic inputs, priors, reference tables, coordinates, or expected invariants.",
    "Replace the shared contract-only runner with a task-specific runner that computes outputs from the raw toy inputs using deterministic local operations and no network access.",
    "Implement 2 to 4 method-shaped computations that finish in under 5 seconds on tiny synthetic inputs and emit machine-readable intermediate QC, not only final deliverables.",
    "Tests must assert numeric or structural invariants of the computation, not only file existence.",
    "Add a short `Starter scope` note to `SKILL.md` stating what is truly computed, what is approximated, and what remains a surrogate for the full public method.",
]


FAMILY_STARTER_GUIDANCE = {
    "spatial_foundation": [
        "Use raw toy inputs that look like reference expression centroids, spot or cell expression, coordinates, and a compact marker table.",
        "Starter operations should include gene intersection, library-size normalize plus log1p, reference-centroid cosine similarity, top-label margin, 3-NN niche majority vote, and marker consistency checks.",
        "Keep the toy problem very small, such as 6 spots, 2 reference cell types, 8 genes, and 3 marker genes.",
        "Do not pre-bake `mapped_labels.tsv`; compute labels and confidence margins from the similarity matrix.",
    ],
    "communication": [
        "Use raw toy inputs that look like cell or group expression matrices, group labels, a ligand-receptor catalog, and an optional adjacency matrix.",
        "Starter operations should include group means, ligand-receptor scoring with `min(ligand_expr, receptor_expr)` or a geometric mean, a small shuffle-based null ranking, and optional adjacency support.",
        "Keep the toy problem very small, such as 8 cells, 3 cell groups, 4 ligand-receptor pairs, and one known positive pair.",
        "For CellPhoneDB or LIANA style tasks, do not pretend to run the full heavy stack; at minimum compute expression filtering, ranking, and a tiny null model.",
    ],
    "multimodal_modeling": [
        "Do not require true totalVI or MultiVI training; implement a lightweight surrogate starter that still computes something real.",
        "Starter operations should normalize RNA and protein separately, concatenate or use a tiny CCA/PCA/SVD latent, build a kNN graph, run label transfer from reference half to query half, and emit batch-mixing plus label-consistency QC.",
        "Keep the toy problem very small, such as 8 cells, 4 RNA genes, 2 proteins, and 2 batches.",
        "The latent representation and QC must be computed from raw toy inputs, not stored directly in a pseudo `.h5mu`.",
    ],
    "trajectory_modeling": [
        "Scope the starter to graph pseudotime plus a tiny velocity surrogate plus a compact fate summary, not a full RNA-velocity stack.",
        "Starter operations should build a kNN graph from toy coordinates or latent features, compute root-to-cell shortest-path or diffusion-style pseudotime, use neighbor differences as a simplified velocity direction, and summarize fate with branch votes or a tiny absorbing Markov chain.",
        "Keep the toy problem very small, such as 10 cells with one trunk, one side branch, one root, and two terminal states.",
        "Tests should check pseudotime monotonicity and higher terminal-state probability on the correct terminal cells.",
    ],
    "tangram_extensions": [
        "Center the starter on gene intersection, mapping weights, and hold-out QC instead of only repeating Tangram prose.",
        "Starter operations should compute cosine-similarity or NNLS-based mapping weights from reference centroids to spatial spots, then use those weights for held-out gene prediction or validation metrics such as marker recovery and assignment entropy.",
        "Keep the toy problem very small, such as 3 spots, 2 cell types, 6 shared genes, and 2 held-out genes.",
        "Outputs must come from computed weights, not from a prewritten `mapping_matrix` inside the toy input.",
    ],
    "spatial_deconvolution": [
        "Define the starter as NNLS-style abundance estimation plus a simple baseline comparison, not as a fake full cell2location or DestVI run.",
        "Starter operations should infer spot-by-cell-type abundances from a reference signature matrix and spot counts, normalize each spot to sum near 1, and compare against a uniform or mean-profile baseline using RMSE or correlation.",
        "Keep the toy problem very small, such as 4 spots, 3 cell types, and 5 genes.",
        "Tests should check per-spot abundance normalization and dominant-type recovery on the target spot.",
    ],
    "bioinformatics_utilities": [
        "Do not hide multiple unrelated tasks behind one abstract template; each task should get a task-specific micro-starter.",
        "`database-query`: use raw toy inputs like alias tables, mock API JSON payloads, and query lists; compute identifier normalization, source merge, and conflict resolution.",
        "`sequence-analysis`: use raw toy inputs like short FASTA records, feature intervals, and variant positions; compute exact or approximate lookup, reverse complement, and a simple primer scan with GC/Tm/product-size rules.",
        "`panel-design`: use raw toy inputs like marker score tables and cell-type coverage matrices; compute marker ranking, redundancy pruning, and coverage balancing.",
    ],
    "spatial_structure": [
        "Use raw toy inputs such as coordinates, expression summaries, neighbor edges, and cluster hints instead of prewritten domain labels.",
        "Starter operations should build the spatial graph, compute compact neighborhood statistics or graph-smoothed scores, assign domains from the resulting structure, and derive a small marker table or enrichment summary.",
        "Keep the toy problem very small, such as 6 to 10 cells or spots across 2 to 3 candidate domains.",
        "Tests should assert that domain assignments are derived from graph structure or statistics, not copied from the toy input.",
    ],
}


COMMON_SPATIAL_CELLS = [
    {"cell_id": "cell_001", "sample_id": "slide_A", "x_coord": 10.0, "y_coord": 5.0},
    {"cell_id": "cell_002", "sample_id": "slide_A", "x_coord": 20.0, "y_coord": 15.0},
    {"cell_id": "cell_003", "sample_id": "slide_A", "x_coord": 30.0, "y_coord": 25.0},
]


def load_manifest() -> dict:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    changed = False
    for task in manifest["tasks"]:
        target_dir = task["target_dir"].replace("experiments/ref_skills/", "experiments/sc_skills/")
        if target_dir != task["target_dir"]:
            task["target_dir"] = target_dir
            changed = True
    for group_name, additions in RESOURCE_GROUP_ADDITIONS.items():
        urls = {item["url"] for item in manifest["resource_groups"].get(group_name, [])}
        for resource in additions:
            if resource["url"] not in urls:
                manifest["resource_groups"].setdefault(group_name, []).append(resource)
                urls.add(resource["url"])
                changed = True
    if changed:
        MANIFEST_PATH.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def task_by_slug(manifest: dict) -> dict[str, dict]:
    return {task["task_slug"]: task for task in manifest["tasks"]}


def select_tasks(manifest: dict, *, task_slugs: list[str], select_all: bool, limit: int | None) -> list[dict]:
    tasks = manifest["tasks"]
    if task_slugs:
        lookup = task_by_slug(manifest)
        selected = [lookup[slug] for slug in task_slugs]
    elif select_all:
        selected = list(tasks)
    else:
        selected = list(tasks[:1])
    if limit is not None:
        selected = selected[:limit]
    return selected


def merged_resources(manifest: dict, task: dict) -> list[dict]:
    merged: list[dict] = []
    seen: set[str] = set()
    for group_name in task["resource_groups"]:
        for resource in manifest["resource_groups"][group_name]:
            if resource["url"] in seen:
                continue
            seen.add(resource["url"])
            merged.append(resource)
    return merged


def skill_front_matter(task: dict) -> str:
    description = (
        f"Use this experiment skill to {task['goal'].rstrip('.')} with public tools and explicit file contracts "
        "that work for Codex, Claude Code, and similar shell-capable agents."
    )
    return "\n".join(
        [
            "---",
            f"name: {task['task_slug']}-portable-skill",
            f"description: {description}",
            "---",
            "",
        ]
    )


def render_skill_md(manifest: dict, task: dict, script_name: str, test_name: str) -> str:
    resources = merged_resources(manifest, task)
    workflow = FAMILY_WORKFLOW[task["family"]]
    deliverables = "\n".join(f"- `{item}`" for item in task["deliverables"])
    methods = "\n".join(f"- {item}" for item in task["default_methods"])
    resources_block = "\n".join(f"- `{item['label']}`" for item in resources)
    focus_terms = ", ".join(task["focus_terms"])
    workflow_block = "\n".join(f"{index}. {step}" for index, step in enumerate(workflow, start=1))
    return (
        skill_front_matter(task)
        + "\n".join(
            [
                "## Purpose",
                "",
                task["goal"],
                "",
                "## Source adaptation",
                "",
                "This experiment skill is derived from `source_reference.md`, but it removes Spatial Agent private tools, hidden wrappers, and notebook-only orchestration.",
                "",
                "## Agent compatibility",
                "",
                "- Compatible with Codex, Claude Code, and similar shell-capable agents.",
                "- Prefer public package CLIs, Python APIs, and explicit file contracts.",
                "- Keep every toy or pseudo-data run reproducible from local files.",
                "",
                "## Focus terms",
                "",
                f"`{focus_terms}`",
                "",
                "## Default methods",
                "",
                methods,
                "",
                "## Recommended resource stack",
                "",
                resources_block,
                "",
                "## Core workflow",
                "",
                workflow_block,
                "",
                "## Deliverables",
                "",
                deliverables,
                "",
                "## Package scaffold",
                "",
                f"- `metadata.yaml` defines the toy contract and required deliverables for `{task['task_slug']}`.",
                f"- `examples/toy_input.json` provides pseudo-data for a deterministic local run.",
                f"- `scripts/{script_name}` is the runnable toy contract entry point.",
                f"- `tests/{test_name}` checks that the local runner writes the declared outputs.",
                "- `assets/README.md` explains what the pseudo-data validates and what it does not.",
                "",
                "## Toy validation",
                "",
                f"- Run `python3 scripts/{script_name} --input examples/toy_input.json --outdir scratch/{task['task_slug']}_toy_run` from this skill directory or with an absolute `--outdir`.",
                f"- Run `python3 -m unittest tests/{test_name}` for the skill-local toy contract check.",
                "- Treat these outputs as toy or pseudo-data validation of the file contract, not as biological evidence.",
                "",
                "## Real-run expectations",
                "",
                "- Replace `examples/toy_input.json` with project or public inputs that preserve the same deliverable interface.",
                "- Pin the exact atlas, database snapshot, model version, and key CLI parameters in the final report or provenance notes.",
                "- Keep stop conditions explicit when feature overlap, spatial metadata, or reference labels are insufficient for a trustworthy run.",
                "",
                "## Starter scope",
                "",
                "- The tiny local starter should compute a small subset of the public method on synthetic or pseudo inputs.",
                "- Keep the boundary explicit between truly computed quantities, lightweight approximations, and placeholders that remain outside starter scope.",
                "- Do not present the starter outputs as equivalent to a full biological analysis on real data.",
                "",
                "## Minimum validation gates",
                "",
                "- All declared deliverables exist and remain machine-readable.",
                "- TSV outputs keep the required columns in `metadata.yaml`.",
                "- Markdown reports contain the named sections used for downstream review.",
                "- `.h5ad` and pseudo-`.h5mu` outputs open successfully and are non-empty.",
                "",
                "## Failure handling",
                "",
                "- Stop early and document blockers when core metadata, coordinates, identifiers, or feature overlap are missing.",
                "- Treat toy outputs as contract validation, not biological evidence.",
                "- Preserve unresolved ambiguity rather than inventing certainty.",
                "",
                "## Hand-off contract",
                "",
                "Before finishing, leave behind the declared deliverables, a concise run summary, and enough provenance for another shell-capable agent to rerun the experiment.",
                "",
            ]
        )
    )


def render_refs_md(manifest: dict, task: dict) -> str:
    lines = ["# References", ""]
    for resource in merged_resources(manifest, task):
        lines.extend(
            [
                f"- {resource['label']}",
                f"  URL: {resource['url']}",
                f"  Note: {resource['note']}",
                "",
            ]
        )
    return "\n".join(lines).rstrip() + "\n"


def render_design_prompt(task: dict) -> str:
    deliverable_block = "\n".join(f"- {item}" for item in task["deliverables"])
    methods = "\n".join(f"- {item}" for item in task["default_methods"])
    hard_rules = "\n".join(f"- {item}" for item in COMMON_PROMPT_RULES)
    family_guidance = "\n".join(f"- {item}" for item in FAMILY_STARTER_GUIDANCE.get(task["family"], []))
    return "\n".join(
        [
            f"Design or refine an experiment-only portable skill for the task: {task['title']}.",
            "",
            "Task goal:",
            task["goal"],
            "",
            "Constraints:",
            f"- Only create or modify files under `{task['target_dir']}`.",
            "- Do not edit `registry/`, `skills/`, `site/`, `README.md`, `experiments.md`, or planning files.",
            "- Keep the result compatible with Codex, Claude Code, and similar shell-capable agents.",
            "- Use `source_reference.md` as decomposition guidance only; remove private Spatial Agent tool dependencies.",
            "- Keep `metadata.yaml`, `examples/`, `scripts/`, `tests/`, and `assets/` aligned with the experiment contract.",
            "- Upgrade the package from a static toy contract toward a small runnable starter that still runs quickly on tiny synthetic or pseudo inputs.",
            "- Prefer method-shaped computations, heuristics, rankings, graph operations, or QC logic over directly dumping fixed files whenever that is feasible on toy data.",
            "- Preserve the declared deliverable names and schemas so repository-level validation still passes.",
            "- If a real public dependency is too heavy for the tiny starter, replace it with a documented lightweight approximation and state the boundary clearly in `SKILL.md` and `assets/README.md`.",
            "",
            "Hard requirements for this refine pass:",
            hard_rules,
            "",
            "Preferred methods:",
            methods,
            "",
            "Family-specific starter guidance:",
            family_guidance,
            "",
            "Required deliverables:",
            deliverable_block,
            "",
            "Refresh the experiment package to keep its docs, scripts, examples, and tests aligned.",
            "Add or refine a tiny runnable example so `scripts/run_exercise.py` performs a meaningful method-shaped starter path instead of only echoing precomputed outputs.",
            "Update `refs.md` if you need additional official docs, papers, notebooks, or GitHub references to justify the starter design.",
            "",
        ]
    )


def make_tsv(columns: list[str], rows: list[dict]) -> dict:
    return {"kind": "tsv", "columns": columns, "rows": rows}


def make_md(title: str, sections: list[tuple[str, list[str]]]) -> dict:
    return {
        "kind": "md",
        "title": title,
        "sections": [{"name": name, "bullets": bullets} for name, bullets in sections],
    }


def make_h5ad(obs_names: list[str], var_names: list[str], matrix: list[list[float]], *, obs: dict | None = None, var: dict | None = None, obsm: dict | None = None) -> dict:
    return {
        "kind": "h5ad",
        "obs_names": obs_names,
        "var_names": var_names,
        "matrix": matrix,
        "obs": obs or {},
        "var": var or {},
        "obsm": obsm or {},
    }


def make_h5mu(modalities: dict[str, dict], *, obsm: dict | None = None) -> dict:
    return {"kind": "h5mu", "modalities": modalities, "obsm": obsm or {}}


def task_example_payloads() -> dict[str, dict]:
    return {
        "annotation": {
            "annotation_table.tsv": make_tsv(
                ["cell_id", "sample_id", "x_coord", "y_coord", "predicted_label", "label_source", "label_confidence", "broad_label", "marker_status", "niche_label", "niche_evidence", "notes"],
                [
                    {**COMMON_SPATIAL_CELLS[0], "predicted_label": "T cell", "label_source": "CellTypist+markers", "label_confidence": 0.94, "broad_label": "immune", "marker_status": "supported", "niche_label": "immune_border", "niche_evidence": "neighbor enrichment", "notes": "pseudo validation"},
                    {**COMMON_SPATIAL_CELLS[1], "predicted_label": "Epithelial", "label_source": "Tangram+markers", "label_confidence": 0.89, "broad_label": "epithelial", "marker_status": "supported", "niche_label": "tumor_core", "niche_evidence": "local composition", "notes": "pseudo validation"},
                ],
            ),
            "marker_evidence.md": make_md("Marker Evidence", [("Run context", ["Pseudo Xenium-like dataset with 3 cells."]), ("Reference source", ["Matched pseudo atlas."]), ("Marker review", ["T cell -> CD3D/CD3E supported.", "Epithelial -> EPCAM/KRT19 supported."]), ("Conflicts and resolutions", ["No major conflicts in toy data."]), ("Remaining uncertainty", ["Niche labels are illustrative only."])]),
            "niche_summary.md": make_md("Niche Summary", [("Run context", ["Pseudo neighborhood graph on three cells."]), ("Spatial graph", ["Used Euclidean neighbor graph from toy coordinates."]), ("Neighborhood validation", ["Immune and epithelial cells are adjacent in the toy layout."]), ("Niche labels", ["Assigned immune_border and tumor_core."]), ("Caveats", ["Toy graph only validates output contract."])]),
        },
        "tissue-niche-annotation": {
            "niche_labels.tsv": make_tsv(
                ["cell_id", "sample_id", "cell_type", "niche_label", "niche_confidence", "dominant_neighbor_type", "boundary_score", "notes"],
                [
                    {"cell_id": "cell_001", "sample_id": "slide_A", "cell_type": "T cell", "niche_label": "immune_aggregate", "niche_confidence": 0.93, "dominant_neighbor_type": "T cell", "boundary_score": 0.08, "notes": "cohesive immune neighborhood"},
                    {"cell_id": "cell_002", "sample_id": "slide_A", "cell_type": "Macrophage", "niche_label": "immune_aggregate", "niche_confidence": 0.84, "dominant_neighbor_type": "T cell", "boundary_score": 0.22, "notes": "myeloid cell embedded in immune neighborhood"},
                    {"cell_id": "cell_003", "sample_id": "slide_A", "cell_type": "Epithelial", "niche_label": "immune_epithelial_interface", "niche_confidence": 0.61, "dominant_neighbor_type": "T cell", "boundary_score": 0.71, "notes": "mixed neighbors near epithelial border"},
                ],
            ),
            "niche_markers.tsv": make_tsv(
                ["niche_label", "dominant_cell_types", "marker_context", "support_score"],
                [
                    {"niche_label": "immune_aggregate", "dominant_cell_types": "T cell;Macrophage", "marker_context": "PTPRC;CXCL9", "support_score": 0.89},
                    {"niche_label": "immune_epithelial_interface", "dominant_cell_types": "Epithelial;T cell", "marker_context": "EPCAM;CXCL10", "support_score": 0.68},
                ],
            ),
            "tissue_niche_report.md": make_md("Tissue Niche Report", [("Run context", ["Pseudo single-cell-resolution spatial niche annotation run."]), ("Spatial graph", ["Neighborhood composition and boundary scores computed from a toy cell graph."]), ("Niche definitions", ["Immune aggregate and immune-epithelial interface were separated from local composition rather than copied from cell types."]), ("Caveats", ["Toy niche labels validate the contract only and are not biological claims."])]),
        },
        "cell-cell-communication": {
            "communication_results.tsv": make_tsv(["source_group", "target_group", "ligand", "receptor", "score", "method", "spatial_support"], [{"source_group": "T_cell", "target_group": "myeloid", "ligand": "CXCL8", "receptor": "CXCR1", "score": 0.91, "method": "LIANA", "spatial_support": "adjacent"}, {"source_group": "myeloid", "target_group": "epithelial", "ligand": "TGFB1", "receptor": "TGFBR2", "score": 0.84, "method": "CellPhoneDB", "spatial_support": "not_tested"}]),
            "priority_pairs.tsv": make_tsv(["rank", "ligand", "receptor", "source_group", "target_group", "priority_reason"], [{"rank": 1, "ligand": "CXCL8", "receptor": "CXCR1", "source_group": "T_cell", "target_group": "myeloid", "priority_reason": "high score plus adjacency"}, {"rank": 2, "ligand": "TGFB1", "receptor": "TGFBR2", "source_group": "myeloid", "target_group": "epithelial", "priority_reason": "cross-method support"}]),
            "interpretation_report.md": make_md("Communication Interpretation", [("Run context", ["Pseudo cell-group communication analysis."]), ("Methods", ["Combined LIANA-like and CellPhoneDB-like toy outputs."]), ("Prioritized interactions", ["CXCL8-CXCR1 ranked highest."]), ("Spatial support", ["Only the top pair had adjacency support in toy data."]), ("Caveats", ["No real expression statistics were computed."])]),
        },
        "cell-deconvolution": {
            "mapping_matrix.h5ad": make_h5ad(["spot_001", "spot_002"], ["T_cell", "myeloid"], [[0.72, 0.28], [0.15, 0.85]], obs={"sample_id": ["slide_A", "slide_A"]}, var={"cell_type": ["T cell", "myeloid"]}, obsm={"spatial": [[5.0, 5.0], [10.0, 10.0]]}),
            "cell_assignment.tsv": make_tsv(["spot_id", "dominant_label", "dominant_score", "secondary_label", "secondary_score"], [{"spot_id": "spot_001", "dominant_label": "T_cell", "dominant_score": 0.72, "secondary_label": "myeloid", "secondary_score": 0.28}, {"spot_id": "spot_002", "dominant_label": "myeloid", "dominant_score": 0.85, "secondary_label": "T_cell", "secondary_score": 0.15}]),
            "deconvolution_report.md": make_md("Deconvolution Report", [("Run context", ["Pseudo Tangram-style constrained mapping."]), ("Input QC", ["Two toy spots and two reference labels."]), ("Assignments", ["Each spot received a dominant label."]), ("Validation", ["Scores sum to one per spot in the toy matrix."]), ("Caveats", ["Pseudo assignment matrix only."])]),
        },
        "cellphonedb-analysis": {
            "cellphonedb_significant_means.tsv": make_tsv(["interacting_pair", "source_group", "target_group", "mean_expression", "pvalue"], [{"interacting_pair": "CXCL8_CXCR1", "source_group": "T_cell", "target_group": "myeloid", "mean_expression": 1.8, "pvalue": 0.01}, {"interacting_pair": "TGFB1_TGFBR2", "source_group": "myeloid", "target_group": "epithelial", "mean_expression": 1.2, "pvalue": 0.03}]),
            "interaction_ranked.tsv": make_tsv(["rank", "interacting_pair", "evidence", "supporting_group_pair"], [{"rank": 1, "interacting_pair": "CXCL8_CXCR1", "evidence": "significant mean + pseudo adjacency", "supporting_group_pair": "T_cell->myeloid"}, {"rank": 2, "interacting_pair": "TGFB1_TGFBR2", "evidence": "significant mean", "supporting_group_pair": "myeloid->epithelial"}]),
            "cellphonedb_report.md": make_md("CellPhoneDB Report", [("Run context", ["Pseudo CellPhoneDB statistical-mode run."]), ("Significant interactions", ["CXCL8_CXCR1 and TGFB1_TGFBR2 retained."]), ("Ranking logic", ["Top rank favors lower p-value plus adjacency support."]), ("Caveats", ["Pseudo statistics only."])]),
        },
        "database-query": {
            "query_results.tsv": make_tsv(["query", "normalized_id", "source", "match_label", "evidence"], [{"query": "TP53", "normalized_id": "ENSG00000141510", "source": "Ensembl", "match_label": "tumor protein p53", "evidence": "canonical gene record"}, {"query": "P04637", "normalized_id": "P04637", "source": "UniProt", "match_label": "Cellular tumor antigen p53", "evidence": "reviewed protein record"}]),
            "source_provenance.md": make_md("Source Provenance", [("Run context", ["Pseudo database lookup workflow."]), ("APIs consulted", ["MyGene.info, Ensembl REST, UniProt."]), ("Normalization", ["TP53 alias normalized before lookup."]), ("Caveats", ["Toy records are deterministic stand-ins."])]),
            "resolution_notes.md": make_md("Resolution Notes", [("Resolved identifiers", ["TP53 -> ENSG00000141510", "P04637 retained as UniProt accession."]), ("Ambiguity handling", ["No ambiguity in toy input."]), ("Follow-up", ["Extend to ClinVar when variant strings are supplied."])]),
        },
        "gene-imputation": {
            "imputed_expression.h5ad": make_h5ad(["cell_001", "cell_002"], ["CXCL13", "MKI67"], [[1.2, 0.1], [0.3, 1.4]], obs={"sample_id": ["slide_A", "slide_A"]}, var={"feature_type": ["gene", "gene"]}, obsm={"spatial": [[0.0, 0.0], [1.0, 1.0]]}),
            "heldout_metrics.tsv": make_tsv(["metric", "value", "higher_is_better"], [{"metric": "pearson_r", "value": 0.82, "higher_is_better": True}, {"metric": "rmse", "value": 0.21, "higher_is_better": False}]),
            "imputation_report.md": make_md("Gene Imputation Report", [("Run context", ["Pseudo Tangram-style gene projection."]), ("Held-out validation", ["Toy Pearson correlation and RMSE recorded."]), ("Interpretation", ["CXCL13 recovered in cell_001."]), ("Caveats", ["Pseudo imputed matrix only."])]),
        },
        "liana-analysis": {
            "liana_rankings.tsv": make_tsv(["rank", "source_group", "target_group", "ligand", "receptor", "aggregate_score"], [{"rank": 1, "source_group": "T_cell", "target_group": "myeloid", "ligand": "CXCL8", "receptor": "CXCR1", "aggregate_score": 0.93}, {"rank": 2, "source_group": "fibroblast", "target_group": "epithelial", "ligand": "COL1A1", "receptor": "ITGB1", "aggregate_score": 0.81}]),
            "resource_overlap.tsv": make_tsv(["pair", "resource_count", "resources"], [{"pair": "CXCL8_CXCR1", "resource_count": 2, "resources": "OmniPath;CellPhoneDB"}, {"pair": "COL1A1_ITGB1", "resource_count": 2, "resources": "OmniPath;LIANA"}]),
            "liana_report.md": make_md("LIANA Report", [("Run context", ["Pseudo LIANA aggregation."]), ("Rankings", ["Two prioritized pairs emitted."]), ("Resource overlap", ["Both top pairs appear in multiple resources."]), ("Caveats", ["No real permutation statistics in toy data."])]),
        },
        "ligand-receptor-discovery": {
            "candidate_pairs.tsv": make_tsv(["candidate_pair", "source_group", "target_group", "support_score"], [{"candidate_pair": "CXCL8_CXCR1", "source_group": "T_cell", "target_group": "myeloid", "support_score": 0.88}, {"candidate_pair": "TGFB1_TGFBR2", "source_group": "myeloid", "target_group": "epithelial", "support_score": 0.79}]),
            "evidence_table.tsv": make_tsv(["candidate_pair", "evidence_type", "evidence_value"], [{"candidate_pair": "CXCL8_CXCR1", "evidence_type": "spatial_support", "evidence_value": "adjacent"}, {"candidate_pair": "TGFB1_TGFBR2", "evidence_type": "resource_support", "evidence_value": "OmniPath"}]),
            "discovery_summary.md": make_md("Ligand-Receptor Discovery Summary", [("Run context", ["Pseudo discovery workflow over candidate pairs."]), ("Prioritized candidates", ["CXCL8_CXCR1 ranked highest."]), ("Evidence summary", ["Spatial and resource support carried through to the evidence table."]), ("Caveats", ["Pseudo discovery only."])]),
        },
        "mapping-validation": {
            "validation_metrics.tsv": make_tsv(["metric", "scope", "value"], [{"metric": "heldout_gene_corr", "scope": "global", "value": 0.77}, {"metric": "marker_recovery", "scope": "label", "value": 0.83}]),
            "holdout_summary.md": make_md("Holdout Summary", [("Run context", ["Pseudo mapping-validation run."]), ("Held-out genes", ["Toy held-out gene correlation recorded."]), ("Label sanity checks", ["Marker recovery above toy threshold."]), ("Caveats", ["No real resampling."])]),
            "mapping_validation_report.md": make_md("Mapping Validation Report", [("Run context", ["Pseudo validation over a precomputed mapping."]), ("Validation metrics", ["Two QC metrics retained."]), ("Interpretation", ["Toy mapping passes the pseudo floor."]), ("Caveats", ["Contract validation only."])]),
        },
        "multimodal-integration": {
            "integrated_latent.h5mu": make_h5mu(
                {
                    "rna": {"obs_names": ["cell_001", "cell_002"], "var_names": ["MS4A1", "CD3D"], "matrix": [[1.0, 0.0], [0.0, 1.0]]},
                    "protein": {"obs_names": ["cell_001", "cell_002"], "var_names": ["CD20", "CD3"], "matrix": [[0.9, 0.1], [0.1, 0.8]]},
                },
                obsm={"X_latent": [[0.1, 0.2], [0.8, 0.7]]},
            ),
            "modality_qc.tsv": make_tsv(["metric", "modality", "value"], [{"metric": "batch_mixing", "modality": "joint", "value": 0.74}, {"metric": "label_consistency", "modality": "joint", "value": 0.88}]),
            "integration_report.md": make_md("Integration Report", [("Run context", ["Pseudo multimodal integration run."]), ("Modalities", ["RNA and protein toy modalities aligned on two cells."]), ("Latent representation", ["Two-dimensional latent embedding stored in pseudo h5mu."]), ("Caveats", ["Pseudo MuData contract only."])]),
        },
        "panel-design": {
            "panel_candidates.tsv": make_tsv(["rank", "target_gene", "panel_role", "rationale"], [{"rank": 1, "target_gene": "EPCAM", "panel_role": "epithelial anchor", "rationale": "high specificity in pseudo atlas"}, {"rank": 2, "target_gene": "PTPRC", "panel_role": "immune anchor", "rationale": "broad immune coverage"}]),
            "panel_rationale.md": make_md("Panel Rationale", [("Run context", ["Pseudo panel-design workflow."]), ("Coverage strategy", ["Balanced epithelial and immune anchors."]), ("Trade-offs", ["Toy panel favors broad classes over fine subtypes."]), ("Caveats", ["No platform chemistry optimization."])]),
            "platform_notes.md": make_md("Platform Notes", [("Platform assumptions", ["Toy notes assume a moderate-plex spatial panel."]), ("Probe constraints", ["Probe feasibility not modeled in pseudo data."]), ("Follow-up", ["Refine with platform-specific constraints in a real run."])]),
        },
        "sequence-analysis": {
            "sequence_summary.tsv": make_tsv(["query_id", "matched_symbol", "matched_accession", "matched_species", "variant_label"], [{"query_id": "query_001", "matched_symbol": "BRAF", "matched_accession": "NM_004333.6", "matched_species": "human", "variant_label": "p.V600E"}]),
            "feature_annotations.tsv": make_tsv(["query_id", "feature_source", "feature_type", "feature_id", "start", "end", "label"], [{"query_id": "query_001", "feature_source": "Ensembl", "feature_type": "exon", "feature_id": "ENSE00003612132", "start": 100, "end": 220, "label": "exon 15"}, {"query_id": "query_001", "feature_source": "UniProt", "feature_type": "domain", "feature_id": "IPR001245", "start": 450, "end": 720, "label": "kinase domain"}]),
            "primer_candidates.tsv": make_tsv(["query_id", "pair_rank", "left_primer", "right_primer", "product_size", "amplicon_target_label"], [{"query_id": "query_001", "pair_rank": 1, "left_primer": "ACGTACGTACGT", "right_primer": "TGCATGCATGCA", "product_size": 180, "amplicon_target_label": "BRAF_V600_region"}]),
        },
        "spatial-deconvolution": {
            "cell_type_abundance.h5ad": make_h5ad(["spot_001", "spot_002"], ["T_cell", "fibroblast"], [[0.6, 0.4], [0.2, 0.8]], obs={"sample_id": ["slide_A", "slide_A"]}, obsm={"spatial": [[1.0, 1.0], [2.0, 2.0]]}),
            "model_comparison.tsv": make_tsv(["model", "metric", "value"], [{"model": "cell2location", "metric": "corr_to_markers", "value": 0.84}, {"model": "DestVI_baseline", "metric": "corr_to_markers", "value": 0.78}]),
            "deconvolution_summary.md": make_md("Deconvolution Summary", [("Run context", ["Pseudo spot-level deconvolution."]), ("Abundance matrix", ["Toy abundance matrix saved as h5ad."]), ("Model comparison", ["Two pseudo models compared on one QC metric."]), ("Caveats", ["No real probabilistic inference."])]),
        },
        "spatial-domain-detection": {
            "domain_labels.tsv": make_tsv(["cell_id", "domain_label", "domain_score"], [{"cell_id": "cell_001", "domain_label": "edge", "domain_score": 0.81}, {"cell_id": "cell_002", "domain_label": "core", "domain_score": 0.88}]),
            "domain_markers.tsv": make_tsv(["domain_label", "marker_gene", "effect_size"], [{"domain_label": "edge", "marker_gene": "CXCL13", "effect_size": 1.2}, {"domain_label": "core", "marker_gene": "KRT19", "effect_size": 1.4}]),
            "domain_detection_report.md": make_md("Domain Detection Report", [("Run context", ["Pseudo spatial-domain detection run."]), ("Domain labels", ["Two domains assigned in toy data."]), ("Markers", ["One illustrative marker per domain."]), ("Caveats", ["No graph neural network training in pseudo mode."])]),
        },
        "spatial-mapping": {
            "mapped_labels.tsv": make_tsv(["observation_id", "x_coord", "y_coord", "primary_label", "primary_score", "score_margin", "secondary_label", "mapping_method", "gene_overlap_count", "gene_overlap_fraction", "marker_review_status", "review_notes"], [{"observation_id": "spot_001", "x_coord": 5.0, "y_coord": 5.0, "primary_label": "T_cell", "primary_score": 0.91, "score_margin": 0.34, "secondary_label": "myeloid", "mapping_method": "Tangram", "gene_overlap_count": 1200, "gene_overlap_fraction": 0.62, "marker_review_status": "supported", "review_notes": "toy mapping"}, {"observation_id": "spot_002", "x_coord": 10.0, "y_coord": 12.0, "primary_label": "epithelial", "primary_score": 0.88, "score_margin": 0.29, "secondary_label": "fibroblast", "mapping_method": "Tangram", "gene_overlap_count": 1200, "gene_overlap_fraction": 0.62, "marker_review_status": "supported", "review_notes": "toy mapping"}]),
            "mapping_scores.tsv": make_tsv(["metric_category", "metric_name", "scope", "target_id", "method", "value", "higher_is_better", "suggested_pass_floor", "qc_status", "notes"], [{"metric_category": "gene_intersection", "metric_name": "shared_genes", "scope": "global", "target_id": "all", "method": "intersection", "value": 1200, "higher_is_better": True, "suggested_pass_floor": 500, "qc_status": "pass", "notes": "toy overlap"}, {"metric_category": "mapping", "metric_name": "score_margin", "scope": "observation", "target_id": "spot_001", "method": "Tangram", "value": 0.34, "higher_is_better": True, "suggested_pass_floor": 0.1, "qc_status": "pass", "notes": "toy margin"}]),
            "spatial_mapping_report.md": make_md("Spatial Mapping Report", [("Run context", ["Pseudo Tangram-style label transfer."]), ("Inputs and preprocessing", ["Reference and spatial toy objects share 1,200 genes."]), ("Gene intersection QC", ["Overlap exceeds the pseudo pass floor."]), ("Tangram mapping", ["Two toy observations were labeled."]), ("Label transfer review", ["Markers support both top labels."]), ("Caveats and next actions", ["Toy mapping only validates the output contract."])]),
        },
        "squidpy-analysis": {
            "spatial_stats.tsv": make_tsv(["statistic", "scope", "value"], [{"statistic": "moran_i", "scope": "gene:CXCL13", "value": 0.41}, {"statistic": "co_occurrence", "scope": "T_cell|myeloid", "value": 0.67}]),
            "neighbor_enrichment.tsv": make_tsv(["label_a", "label_b", "enrichment_score"], [{"label_a": "T_cell", "label_b": "myeloid", "enrichment_score": 1.4}, {"label_a": "epithelial", "label_b": "fibroblast", "enrichment_score": 1.2}]),
            "squidpy_report.md": make_md("Squidpy Report", [("Run context", ["Pseudo Squidpy neighborhood analysis."]), ("Spatial graph", ["Toy Euclidean graph across four cells."]), ("Statistics", ["Two spatial statistics recorded."]), ("Caveats", ["No real permutation testing."])]),
        },
        "trajectory-inference": {
            "trajectory_coordinates.tsv": make_tsv(["cell_id", "pseudotime", "lineage_label"], [{"cell_id": "cell_001", "pseudotime": 0.05, "lineage_label": "root"}, {"cell_id": "cell_002", "pseudotime": 0.55, "lineage_label": "branch_A"}, {"cell_id": "cell_003", "pseudotime": 0.92, "lineage_label": "terminal_A"}]),
            "fate_probabilities.tsv": make_tsv(["cell_id", "terminal_state", "probability"], [{"cell_id": "cell_001", "terminal_state": "terminal_A", "probability": 0.22}, {"cell_id": "cell_002", "terminal_state": "terminal_A", "probability": 0.64}, {"cell_id": "cell_003", "terminal_state": "terminal_A", "probability": 0.97}]),
            "trajectory_report.md": make_md("Trajectory Report", [("Run context", ["Pseudo trajectory inference run."]), ("Coordinates", ["Toy pseudotime values span root to terminal state."]), ("Fate summaries", ["One terminal fate probability table recorded."]), ("Caveats", ["No real velocity or fate model fit."])]),
        },
    }


def metadata_from_example(task: dict, example_payload: dict) -> dict:
    deliverables = []
    for path in task["deliverables"]:
        spec = example_payload[path]
        item = {"path": path, "kind": spec["kind"]}
        if spec["kind"] == "tsv":
            item["required_columns"] = list(spec["columns"])
        elif spec["kind"] == "md":
            item["required_sections"] = [section["name"] for section in spec["sections"]]
        deliverables.append(item)
    return {
        "task_slug": task["task_slug"],
        "title": task["title"],
        "family": task["family"],
        "experiment_only": True,
        "source_reference": task["source_reference"],
        "target_dir": task["target_dir"],
        "focus_terms": task["focus_terms"],
        "default_methods": task["default_methods"],
        "deliverables": deliverables,
    }


def script_name(task_slug: str) -> str:
    return f"run_{task_slug.replace('-', '_')}.py"


def test_name(task_slug: str) -> str:
    return f"test_{task_slug.replace('-', '_')}.py"


def render_runner_wrapper(task_slug: str) -> str:
    return "\n".join(
        [
            "#!/usr/bin/env python3",
            "from __future__ import annotations",
            "",
            "from pathlib import Path",
            "import sys",
            "",
            "ROOT = Path(__file__).resolve().parents[4]",
            "sys.path.insert(0, str(ROOT))",
            "",
            "from experiments.sc_skills.shared.contract_runner import main_for_skill",
            "",
            "",
            "if __name__ == \"__main__\":",
            "    raise SystemExit(main_for_skill(Path(__file__).resolve().parents[1]))",
            "",
        ]
    )


def render_run_exercise_wrapper(task_slug: str) -> str:
    return "\n".join(
        [
            "#!/usr/bin/env python3",
            "from __future__ import annotations",
            "",
            "import argparse",
            "from pathlib import Path",
            "import subprocess",
            "import sys",
            "",
            "SKILL_DIR = Path(__file__).resolve().parents[1]",
            f"SCRIPT = SKILL_DIR / 'scripts' / '{script_name(task_slug)}'",
            "VALIDATOR = SKILL_DIR / 'scripts' / 'validate_outputs.py'",
            "INPUT = SKILL_DIR / 'examples' / 'toy_input.json'",
            "",
            "parser = argparse.ArgumentParser()",
            "parser.add_argument('--outdir', type=Path, default=SKILL_DIR / 'scratch' / 'toy_run')",
            "args = parser.parse_args()",
            "args.outdir.mkdir(parents=True, exist_ok=True)",
            "completed = subprocess.run(",
            "    [sys.executable, str(SCRIPT), '--input', str(INPUT), '--outdir', str(args.outdir)],",
            "    cwd=SKILL_DIR.parents[2],",
            "    check=False,",
            "    text=True,",
            ")",
            "if completed.returncode != 0:",
            "    raise SystemExit(completed.returncode)",
            "validated = subprocess.run(",
            "    [sys.executable, str(VALIDATOR), '--outdir', str(args.outdir)],",
            "    cwd=SKILL_DIR.parents[2],",
            "    check=False,",
            "    text=True,",
            ")",
            "raise SystemExit(validated.returncode)",
            "",
        ]
    )


def render_validate_outputs_wrapper() -> str:
    return "\n".join(
        [
            "#!/usr/bin/env python3",
            "from __future__ import annotations",
            "",
            "import argparse",
            "from pathlib import Path",
            "import sys",
            "",
            "ROOT = Path(__file__).resolve().parents[4]",
            "sys.path.insert(0, str(ROOT))",
            "",
            "from experiments.sc_skills.shared.contract_runner import validate_outputs",
            "",
            "",
            "parser = argparse.ArgumentParser()",
            "parser.add_argument('--outdir', type=Path, required=True)",
            "args = parser.parse_args()",
            "validate_outputs(Path(__file__).resolve().parents[1], args.outdir)",
            "",
        ]
    )


def render_skill_local_test(task_slug: str) -> str:
    return "\n".join(
        [
            "from __future__ import annotations",
            "",
            "import json",
            "import subprocess",
            "import tempfile",
            "import unittest",
            "from pathlib import Path",
            "",
            "ROOT = Path(__file__).resolve().parents[4]",
            "SKILL_DIR = Path(__file__).resolve().parents[1]",
            "SCRIPT = SKILL_DIR / 'scripts' / 'run_exercise.py'",
            "",
            "",
            "class SkillLocalToyRunTests(unittest.TestCase):",
            "    def test_toy_run_writes_declared_outputs(self) -> None:",
            "        with tempfile.TemporaryDirectory() as tmpdir:",
            "            metadata = json.loads((SKILL_DIR / 'metadata.yaml').read_text(encoding='utf-8'))",
            "            completed = subprocess.run(",
            "                ['python3', str(SCRIPT), '--outdir', tmpdir],",
            "                cwd=ROOT,",
            "                check=False,",
            "                capture_output=True,",
            "                text=True,",
            "            )",
            "            self.assertEqual(completed.returncode, 0, completed.stderr)",
            "            for deliverable in metadata['deliverables']:",
            "                self.assertTrue((Path(tmpdir) / deliverable['path']).exists(), deliverable['path'])",
            "            self.assertTrue((Path(tmpdir) / 'run_summary.json').exists())",
            "",
            "",
            "if __name__ == '__main__':",
            "    unittest.main()",
            "",
        ]
    )


def render_contract_test() -> str:
    return "\n".join(
        [
            "from __future__ import annotations",
            "",
            "import json",
            "import unittest",
            "from pathlib import Path",
            "",
            "SKILL_DIR = Path(__file__).resolve().parents[1]",
            "",
            "",
            "class SkillContractMetadataTests(unittest.TestCase):",
            "    def test_metadata_matches_example_deliverables(self) -> None:",
            "        metadata = json.loads((SKILL_DIR / 'metadata.yaml').read_text(encoding='utf-8'))",
            "        example = json.loads((SKILL_DIR / 'examples' / 'toy_input.json').read_text(encoding='utf-8'))",
            "        expected = sorted(item['path'] for item in metadata['deliverables'])",
            "        observed = sorted(example['deliverables'].keys())",
            "        self.assertEqual(expected, observed)",
            "",
            "",
            "if __name__ == '__main__':",
            "    unittest.main()",
            "",
        ]
    )


def render_assets_readme(task: dict) -> str:
    return "\n".join(
        [
            "# Assets",
            "",
            f"This directory belongs to the experiment skill `{task['task_slug']}`.",
            "",
            "The current hardening pass uses deterministic toy or pseudo-data only. No real biological claims should be drawn from these assets.",
            "",
            "## What the toy assets validate",
            "",
            "- The runnable contract between `examples/`, `scripts/`, `metadata.yaml`, and the declared deliverables.",
            "- Basic machine-readability of `.tsv`, `.md`, `.h5ad`, and pseudo-`.h5mu` artifacts.",
            "- Repeatable output structure that a general coding agent can regenerate locally.",
            "",
            "## What the toy assets do not validate",
            "",
            "- Method-specific biological correctness on real spatial or single-cell data.",
            "- External database reachability, model calibration, or atlas selection quality.",
            "- Statistical claims beyond the deterministic pseudo-data contract.",
            "",
            "## Upgrade path",
            "",
            "- Replace the toy input with a pinned real dataset or public benchmark slice.",
            "- Record parameter choices and provenance in the report deliverables and `refs.md`.",
            "- Keep the same deliverable names and schemas when promoting the experiment to a stronger runnable starter.",
            "",
        ]
    )


def render_resource_groups_md(manifest: dict) -> str:
    lines = ["# Resource Groups", "", "Canonical grouped resources for `experiments/sc_skills/`.", ""]
    for group_name, resources in sorted(manifest["resource_groups"].items()):
        lines.append(f"## {group_name}")
        lines.append("")
        for resource in resources:
            lines.append(f"- {resource['label']}")
            lines.append(f"  URL: {resource['url']}")
            lines.append(f"  Note: {resource['note']}")
            lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def render_inventory_md(manifest: dict) -> str:
    lines = [
        "# Task Inventory",
        "",
        "This inventory maps every Spatial Agent reference under [`ref/skill/`](/home/shuaikes/projects/agent/SciSkillUniverse/ref/skill) to an experiment-only portable skill target under [`experiments/sc_skills/`](/home/shuaikes/projects/agent/SciSkillUniverse/experiments/sc_skills).",
        "",
        "| Task slug | Title | Family | Deliverables |",
        "| --- | --- | --- | --- |",
    ]
    for task in manifest["tasks"]:
        lines.append(f"| {task['task_slug']} | {task['title']} | {task['family']} | {', '.join(task['deliverables'])} |")
    lines.extend(
        [
            "",
            f"- The machine-readable source of truth for this inventory is [`batch_design_manifest.json`](/home/shuaikes/projects/agent/SciSkillUniverse/experiments/sc_skills/batch_design_manifest.json).",
            "- Every experiment package is expected to carry docs plus `metadata.yaml`, `scripts/`, `examples/`, `tests/`, and `assets/` after hardening.",
            "",
        ]
    )
    return "\n".join(lines)


def render_readme() -> str:
    return "\n".join(
        [
            "# sc_skills Experiments",
            "",
            "This directory is the experiment-only workspace for portable single-cell and spatial skills derived from the Spatial Agent task references in `ref/skill/`.",
            "",
            "## Scope",
            "",
            "- Source references live in [`ref/skill/`](/home/shuaikes/projects/agent/SciSkillUniverse/ref/skill).",
            "- Generated experiment skills stay under [`experiments/sc_skills/`](/home/shuaikes/projects/agent/SciSkillUniverse/experiments/sc_skills).",
            "- These experiment skills are not production skills and are not registered in `registry/`.",
            "- The batch design flow can optionally route a selected task through the framework `design-skill` command.",
            "",
            "## Files",
            "",
            "- [`task_inventory.md`](/home/shuaikes/projects/agent/SciSkillUniverse/experiments/sc_skills/task_inventory.md)",
            "- [`resource_groups.md`](/home/shuaikes/projects/agent/SciSkillUniverse/experiments/sc_skills/resource_groups.md)",
            "- [`batch_design_manifest.json`](/home/shuaikes/projects/agent/SciSkillUniverse/experiments/sc_skills/batch_design_manifest.json)",
            "",
            "## Batch usage",
            "",
            "```bash",
            "python3 scripts/batch_design_experiment_skills.py --print-only --limit 3",
            "python3 scripts/batch_design_experiment_skills.py --task annotation",
            "python3 scripts/batch_design_experiment_skills.py --all",
            "python3 scripts/batch_design_experiment_skills.py --all --run-framework --verification-mode validate",
            "python3 scripts/validate_sc_skill_experiments.py --all --json-out scratch/reviews/sc_skills_validation.json --markdown-out scratch/reviews/sc_skills_validation.md",
            "python3 scripts/run_ref_skill_framework_batches.py run --print-plan",
            "```",
            "",
            "Each task is expected to include documentation plus runnable toy scaffolding under `scripts/`, `examples/`, `tests/`, and `assets/`.",
            "",
        ]
    )


def materialize_task(manifest: dict, task: dict, example_payloads: dict[str, dict]) -> Path:
    target_dir = ROOT / task["target_dir"]
    target_dir.mkdir(parents=True, exist_ok=True)
    source_path = ROOT / task["source_reference"]
    source_text = source_path.read_text(encoding="utf-8")
    (target_dir / "source_reference.md").write_text(source_text, encoding="utf-8")
    skill_script = script_name(task["task_slug"])
    skill_test = test_name(task["task_slug"])
    (target_dir / "SKILL.md").write_text(render_skill_md(manifest, task, skill_script, skill_test), encoding="utf-8")
    (target_dir / "refs.md").write_text(render_refs_md(manifest, task), encoding="utf-8")
    (target_dir / "design_prompt.md").write_text(render_design_prompt(task), encoding="utf-8")

    example_payload = {"run_label": f"toy-{task['task_slug']}", "deliverables": example_payloads[task["task_slug"]]}
    metadata = metadata_from_example(task, example_payload["deliverables"])

    (target_dir / "metadata.yaml").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (target_dir / "examples").mkdir(exist_ok=True)
    (target_dir / "examples" / "toy_input.json").write_text(json.dumps(example_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (target_dir / "scripts").mkdir(exist_ok=True)
    script_path = target_dir / "scripts" / skill_script
    script_path.write_text(render_runner_wrapper(task["task_slug"]), encoding="utf-8")
    script_path.chmod(0o755)
    exercise_path = target_dir / "scripts" / "run_exercise.py"
    exercise_path.write_text(render_run_exercise_wrapper(task["task_slug"]), encoding="utf-8")
    exercise_path.chmod(0o755)
    validate_path = target_dir / "scripts" / "validate_outputs.py"
    validate_path.write_text(render_validate_outputs_wrapper(), encoding="utf-8")
    validate_path.chmod(0o755)
    (target_dir / "tests").mkdir(exist_ok=True)
    (target_dir / "tests" / skill_test).write_text(render_skill_local_test(task["task_slug"]), encoding="utf-8")
    (target_dir / "tests" / "test_contract.py").write_text(render_contract_test(), encoding="utf-8")
    (target_dir / "assets").mkdir(exist_ok=True)
    (target_dir / "assets" / "README.md").write_text(render_assets_readme(task), encoding="utf-8")
    return target_dir


def sync_root_docs(manifest: dict) -> None:
    (EXPERIMENT_ROOT / "README.md").write_text(render_readme() + "\n", encoding="utf-8")
    (EXPERIMENT_ROOT / "resource_groups.md").write_text(render_resource_groups_md(manifest), encoding="utf-8")
    (EXPERIMENT_ROOT / "task_inventory.md").write_text(render_inventory_md(manifest) + "\n", encoding="utf-8")


def run_framework(
    task: dict,
    *,
    verification_mode: str,
    model: str | None = None,
    reasoning_effort: str | None = None,
    profile: str | None = None,
    codex_bin: str | None = None,
    full_auto: bool = False,
    extra_context: str | None = None,
) -> subprocess.CompletedProcess[str]:
    prompt_path = ROOT / task["target_dir"] / "design_prompt.md"
    extra_context = (
        f"Strictly stay inside {task['target_dir']}. "
        "Do not touch production skills, registries, site assets, or global docs. "
        "This is an experiment-only sc-skill port. "
        "Upgrade the current toy contract into a small runnable starter closer to the cited public methods while keeping the declared deliverable interface stable."
        + (f" {extra_context.strip()}" if extra_context else "")
    )
    command = [
        "python3",
        "scripts/sciskill_framework.py",
    ]
    if codex_bin:
        command.extend(["--codex-bin", codex_bin])
    if model:
        command.extend(["--model", model])
    if reasoning_effort:
        command.extend(["--reasoning-effort", reasoning_effort])
    if profile:
        command.extend(["--profile", profile])
    if full_auto:
        command.append("--full-auto")
    command.extend(
        [
        "design-skill",
        "--prompt",
        prompt_path.read_text(encoding="utf-8"),
        "--verification-mode",
        verification_mode,
        "--label",
        f"sc-skill-{task['task_slug']}",
        "--extra-context",
        extra_context,
        ]
    )
    for focus_term in task["focus_terms"]:
        command.extend(["--focus-term", focus_term])
    return subprocess.run(command, cwd=ROOT, check=False, capture_output=True, text=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task", action="append", default=[], help="Task slug from the batch manifest.")
    parser.add_argument("--all", action="store_true", help="Select all manifest tasks.")
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--print-only", action="store_true", help="Print the selected tasks and exit.")
    parser.add_argument("--run-framework", action="store_true", help="After materialization, call the framework design-skill entry for each selected task.")
    parser.add_argument("--verification-mode", default="none", choices=["none", "validate", "standard", "full", "audit"])
    parser.add_argument("--framework-model", default=None)
    parser.add_argument("--framework-reasoning-effort", default=None)
    parser.add_argument("--framework-profile", default=None)
    parser.add_argument("--framework-codex-bin", default=None)
    parser.add_argument("--framework-full-auto", action="store_true")
    parser.add_argument("--framework-extra-context", default=None)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    manifest = load_manifest()
    sync_root_docs(manifest)
    selected = select_tasks(manifest, task_slugs=args.task, select_all=args.all, limit=args.limit)

    if args.print_only:
        print(json.dumps([task["task_slug"] for task in selected], indent=2))
        return 0

    example_payloads = task_example_payloads()
    materialized: list[str] = []
    framework_runs: list[dict] = []
    for task in selected:
        target_dir = materialize_task(manifest, task, example_payloads)
        materialized.append(str(target_dir.relative_to(ROOT)))
        if args.run_framework:
            completed = run_framework(
                task,
                verification_mode=args.verification_mode,
                model=args.framework_model,
                reasoning_effort=args.framework_reasoning_effort,
                profile=args.framework_profile,
                codex_bin=args.framework_codex_bin,
                full_auto=args.framework_full_auto,
                extra_context=args.framework_extra_context,
            )
            framework_runs.append(
                {
                    "task_slug": task["task_slug"],
                    "returncode": completed.returncode,
                    "stdout_tail": completed.stdout.splitlines()[-20:],
                    "stderr_tail": completed.stderr.splitlines()[-20:],
                }
            )

    payload = {"materialized": materialized, "framework_runs": framework_runs}
    print(json.dumps(payload, indent=2))
    return 0 if all(item.get("returncode", 0) == 0 for item in framework_runs) else 1


if __name__ == "__main__":
    raise SystemExit(main())
