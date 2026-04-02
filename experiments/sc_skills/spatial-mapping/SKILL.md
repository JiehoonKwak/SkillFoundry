---
name: spatial-mapping-portable-skill
description: Use this experiment skill to map scRNA-seq reference cells or labels onto spatial observations with explicit preprocessing, gene intersection, and QC steps using portable local workflows that work for Codex, Claude Code, and similar shell-capable agents.
---
## Purpose

Map scRNA-seq reference cells or labels onto spatial observations with explicit preprocessing, gene intersection, and QC steps.

## Source adaptation

This experiment skill is derived from `source_reference.md`, but it removes Spatial Agent private tools, hidden wrappers, and notebook-only orchestration.

## Agent compatibility

- Compatible with Codex, Claude Code, and similar shell-capable agents.
- Prefer public package CLIs, Python APIs, and explicit file contracts.
- Keep every toy or pseudo-data run reproducible from local files.

## Focus terms

`spatial mapping, Tangram, label transfer`

## Default methods

- Tangram mapping
- gene intersection QC
- label transfer review

## Recommended resource stack

- Tangram documentation for `pp_adatas`, `map_cells_to_space`, and mapping tutorials
- Tangram Nature Methods paper
- Scanpy `normalize_total` and `log1p` docs
- Squidpy `spatial_neighbors` docs
- CellTypist docs for best-match confidence and majority-voting review patterns
- Spatial integration benchmark notes for Tangram positioning

## Core workflow

1. Inspect the reference/spatial gene overlap, coordinate fields, and label surface before attempting transfer.
2. Intersect genes explicitly, record dropped features, and normalize each profile with a `normalize_total` plus `log1p`-shaped preprocessing step.
3. Score each spatial observation against reference centroids with cosine similarity and keep the full similarity matrix in machine-readable QC.
4. Review the top label with score margins, marker consistency, and 3-nearest-neighbor niche voting before finalizing mapped labels.

## Deliverables

- `mapped_labels.tsv`
- `mapping_scores.tsv`
- `spatial_mapping_report.md`

## Package scaffold

- `metadata.yaml` defines the stable deliverable contract plus required QC sidecars for `spatial-mapping`.
- `examples/toy_input.json` contains raw synthetic reference centroids, spatial counts, coordinates, marker tables, and expected invariants.
- `scripts/run_spatial_mapping.py` computes gene intersection, normalize+log1p preprocessing, centroid cosine scores, top-label margins, 3-NN niche votes, and marker review output without network access.
- `scripts/validate_outputs.py` rechecks the declared deliverables plus task-specific numeric invariants.
- `tests/test_contract.py` checks that the toy input remains raw-input-focused and aligned with the metadata contract.
- `tests/test_spatial_mapping.py` asserts structural and numeric invariants for the starter computations.
- `assets/README.md` explains the computed versus approximated starter boundary.

## Toy validation

- Run `python3 scripts/run_spatial_mapping.py --input examples/toy_input.json --outdir scratch/spatial-mapping_toy_run` from this skill directory or with an absolute `--outdir`.
- Run `python3 scripts/validate_outputs.py --outdir scratch/spatial-mapping_toy_run --input examples/toy_input.json` to recheck deliverables and toy invariants.
- Run `python3 -m unittest tests/test_contract.py tests/test_spatial_mapping.py` for the skill-local test pass.
- Inspect `intermediate_qc.json`, `gene_intersection_qc.tsv`, `marker_qc.tsv`, and `niche_qc.tsv` for the computed intermediate checks.

## Real-run expectations

- Replace `examples/toy_input.json` with project or public inputs that preserve the same deliverable interface.
- Pin the exact atlas slice, label vocabulary, normalization choices, and mapping parameters in the final report or provenance notes.
- Keep stop conditions explicit when feature overlap, coordinates, or label evidence are insufficient for a trustworthy run.

## Starter scope

- Computed: gene intersection, dropped-feature QC, library-size normalization with `log1p`, reference-centroid cosine similarity, top-label margins, 3-NN niche majority voting, and marker consistency checks on tiny synthetic inputs.
- Approximated: Tangram is represented by a centroid-level cosine surrogate rather than the full PyTorch optimizer, density priors, or projected-gene workflow.
- Surrogate: the toy reference centroids, marker rules, and review thresholds are synthetic, so the starter validates workflow shape and QC logic rather than real biological performance.

## Minimum validation gates

- All declared deliverables exist and remain machine-readable.
- `mapping_scores.tsv` keeps the required columns in `metadata.yaml`.
- Markdown deliverables contain the named sections used for downstream review.
- `intermediate_qc.json`, `gene_intersection_qc.tsv`, `marker_qc.tsv`, and `niche_qc.tsv` stay consistent with the fitted labels and margins.

## Failure handling

- Stop early and document blockers when core metadata, coordinates, identifiers, or feature overlap are missing.
- Treat toy outputs as contract validation, not biological evidence.
- Preserve unresolved ambiguity rather than inventing certainty.

## Hand-off contract

Before finishing, leave behind the declared deliverables, the intermediate QC files, a concise run summary, and enough provenance for another shell-capable agent to rerun the experiment.
