---
name: spatial-domain-detection-portable-skill
description: Use this experiment skill to detect spatial regions or domains with graph-based models and explain the evidence behind the chosen partition with public tools and explicit file contracts that work for Codex, Claude Code, and similar shell-capable agents.
---
## Purpose

Detect spatial regions or domains with graph-based models and explain the evidence behind the chosen partition.

## Source adaptation

This experiment skill is derived from `source_reference.md`, but it removes Spatial Agent private tools, hidden wrappers, and notebook-only orchestration.

## Agent compatibility

- Compatible with Codex, Claude Code, and similar shell-capable agents.
- Prefer public package CLIs, Python APIs, and explicit file contracts.
- Keep every toy or pseudo-data run reproducible from local files.

## Focus terms

`spatial domains, SpaGCN, GraphST`

## Default methods

- SpaGCN
- GraphST
- domain marker interpretation

## Recommended resource stack

- `SpaGCN paper and repository`
- `GraphST paper, repository, and docs`
- `Squidpy spatial-neighbors documentation`
- `Scanpy rank_genes_groups documentation`

## Core workflow

1. Build compact spatial graphs from explicit coordinates, with a histology-aware variant for the SpaGCN-like path and a coordinate-only variant for the GraphST-like path.
2. Score candidate domains from lightweight domain signatures on normalized toy expression.
3. Compare a graph-smoothed SpaGCN-like partition against a graph-diffused GraphST-like embedding partition.
4. Emit domain labels, marker effects, and machine-readable QC that explain why the selected partition won.

## Deliverables

- `domain_labels.tsv`
- `domain_markers.tsv`
- `domain_detection_report.md`

## Package scaffold

- `metadata.yaml` defines the deliverables and starter QC files.
- `examples/toy_input.json` holds only raw synthetic spots, coordinates, histology proxies, signatures, and expected invariants.
- `scripts/run_spatial_domain_detection.py` computes the deterministic starter path locally.
- `scripts/validate_outputs.py` checks deliverables plus graph/QC invariants.
- `tests/test_spatial_domain_detection.py` asserts graph-derived assignment and marker invariants.
- `assets/README.md` records the precise starter boundary against the full public methods.

## Toy validation

- Run `python3 scripts/run_spatial_domain_detection.py --input examples/toy_input.json --outdir scratch/spatial-domain-detection_toy_run`.
- Run `python3 scripts/validate_outputs.py --outdir scratch/spatial-domain-detection_toy_run`.
- Run `python3 -m unittest tests/test_contract.py tests/test_spatial_domain_detection.py`.
- Treat these outputs as a deterministic starter, not biological evidence.

## Real-run expectations

- Replace `examples/toy_input.json` with project or public inputs that preserve the same deliverable interface.
- Pin the exact dataset slice, graph parameters, and any atlas or histology preprocessing in the final report.
- Keep stop conditions explicit when coordinates, expression summaries, or candidate domain priors are too weak for a trustworthy run.

## Starter scope

- Truly computed: spatial neighbor graphs, graph-smoothed domain scores, graph-diffused embeddings with deterministic clustering, marker effect sizes, and method-comparison QC on tiny toy inputs.
- Approximated: SpaGCN histology integration is reduced to RGB-like proxy affinities, and GraphST self-supervision is reduced to one-step diffusion plus SVD and k-means.
- Surrogate only: no trained GCN, no contrastive learning, no real H&E feature extraction, and no full statistical DE workflow on real spatial transcriptomics data.

## Minimum validation gates

- All declared deliverables exist and remain machine-readable.
- TSV outputs keep the required columns in `metadata.yaml`.
- Markdown reports contain the named sections used for downstream review.
- QC outputs include both graph variants, per-spot method scores, a selected-method summary, and row-normalized graph weights.

## Failure handling

- Stop early and document blockers when coordinates, identifiers, or signature genes are missing.
- Treat toy outputs as contract validation, not biological evidence.
- Preserve unresolved ambiguity rather than inventing certainty.

## Hand-off contract

Before finishing, leave behind the declared deliverables, the QC artifacts, and enough provenance for another shell-capable agent to rerun the starter exactly.
