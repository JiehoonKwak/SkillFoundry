---
name: squidpy-analysis-portable-skill
description: Use this experiment skill to run portable Squidpy-style spatial graphs, neighborhood enrichment, co-occurrence, and Moran's I with lightweight image-aware surrogates and explicit file contracts that work for Codex, Claude Code, and similar shell-capable agents.
---
## Purpose

Run portable Squidpy-style spatial graphs, neighborhood enrichment, and image-aware spatial statistics.

## Source adaptation

This experiment skill is derived from `source_reference.md`, but it removes private tool dependencies and notebook-only orchestration in favor of public, local, deterministic steps.

## Agent compatibility

- Compatible with Codex, Claude Code, and similar shell-capable agents.
- Prefer public Python APIs, explicit file contracts, and toy inputs that another shell-capable agent can rerun locally.
- Keep every starter computation deterministic and fast enough for local validation.

## Focus terms

`squidpy, spatial neighbors, neighborhood enrichment, co-occurrence, Moran's I`

## Default methods

- spatial neighbors
- co-occurrence
- Moran's I

## Recommended resource stack

- `Squidpy graph API docs`
- `Squidpy image feature docs`
- `Squidpy ImageContainer tutorial`
- `Squidpy paper`
- `Squidpy GitHub repository`

## Core workflow

1. Build an image-aware spatial KNN graph from explicit coordinates plus compact per-spot image summary features.
2. Score candidate domains from raw expression and image priors, then smooth the scores over the graph to derive neighborhood labels for downstream statistics.
3. Compute deterministic neighborhood enrichment and co-occurrence summaries from the derived labels.
4. Compute Moran's I across toy gene summaries and toy image features, then write deliverables plus machine-readable QC.

## Deliverables

- `spatial_stats.tsv`
- `neighbor_enrichment.tsv`
- `squidpy_report.md`

## Package scaffold

- `metadata.yaml` defines the stable deliverable contract and starter QC files.
- `examples/toy_input.json` holds only raw synthetic spots, coordinates, counts, image feature summaries, priors, and expected invariants.
- `scripts/run_squidpy_analysis.py` computes the deterministic starter path locally.
- `scripts/validate_outputs.py` checks deliverables plus graph, label-shift, and statistic invariants.
- `tests/test_squidpy_analysis.py` asserts graph-derived label shifts and numeric spatial-statistic invariants.
- `assets/README.md` records the precise starter boundary against full Squidpy workflows.

## Toy validation

- Run `python3 scripts/run_squidpy_analysis.py --input examples/toy_input.json --outdir scratch/squidpy-analysis_toy_run`.
- Run `python3 scripts/validate_outputs.py --outdir scratch/squidpy-analysis_toy_run --input examples/toy_input.json`.
- Run `python3 -m unittest tests/test_contract.py tests/test_squidpy_analysis.py`.
- Treat these outputs as a deterministic starter, not biological evidence.

## Real-run expectations

- Replace `examples/toy_input.json` with project or public inputs that preserve the same deliverable interface.
- Pin the exact graph parameters, image-summary provenance, and feature preprocessing in the final report.
- Stop early when coordinates, identifiers, or image summaries are too weak to support trustworthy spatial statistics.

## Starter scope

- Truly computed: image-aware spatial neighbors, graph-smoothed label assignment from raw toy counts and image summaries, deterministic neighborhood enrichment, co-occurrence within a fixed distance threshold, and Moran's I over toy gene/image features.
- Approximated: the image-aware path uses precomputed per-spot summary features instead of Squidpy `ImageContainer` crops, and neighborhood enrichment is a log2 observed-versus-expected surrogate rather than a permutation z-score.
- Surrogate only: no real AnnData object, no full Squidpy permutation testing, no significance calibration, and no direct histology patch extraction or segmentation on real images.

## Minimum validation gates

- All declared deliverables exist and remain machine-readable.
- TSV outputs keep the required columns in `metadata.yaml`.
- Markdown reports contain the named sections used for downstream review.
- QC outputs include the spatial graph, graph-derived label scores, pairwise enrichment/co-occurrence statistics, and `run_summary.json`.

## Failure handling

- Stop early and document blockers when coordinates, identifiers, raw counts, or image summaries are missing.
- Treat toy outputs as contract validation, not biological evidence.
- Preserve unresolved ambiguity rather than inventing certainty.

## Hand-off contract

Before finishing, leave behind the declared deliverables, the QC artifacts, and enough provenance for another shell-capable agent to rerun the starter exactly.
