---
name: multimodal-integration-portable-skill
description: Use this experiment skill to integrate RNA, protein, and chromatin modalities with explicit model choice, latent-space QC, and label-transfer outputs using a small deterministic starter that works in Codex, Claude Code, and similar shell-capable agents.
---
## Purpose

Integrate RNA, protein, and chromatin modalities with explicit model choice, latent-space QC, and label-transfer outputs.

## Source adaptation

This experiment skill is derived from `source_reference.md`, but it removes private Spatial Agent tools, hidden wrappers, and notebook-only orchestration.

## Agent compatibility

- Compatible with Codex, Claude Code, and similar shell-capable agents.
- Uses only local files plus public Python dependencies already present in the repository environment.
- Keeps the deliverable contract explicit and reproducible from shell commands.

## Focus terms

`multimodal integration, totalVI, MultiVI, latent QC, label transfer`

## Starter methods

- `multivi_style`: weighted RNA + chromatin latent with protein as an auxiliary block.
- `totalvi_style`: weighted RNA + protein latent with chromatin as an auxiliary block.
- Shared downstream steps: weighted SVD latent construction, kNN graph QC, and reference-to-query label transfer.

## Deliverables

- `integrated_latent.h5mu`
- `modality_qc.tsv`
- `integration_report.md`

## Starter scope

- Truly computed: modality-wise normalization, weighted tri-modal latent coordinates, kNN graph structure, batch-mixing QC, label-consistency QC, and kNN label transfer on raw toy counts.
- Approximated: totalVI and MultiVI are represented by deterministic weighted-fusion surrogates rather than trained variational models.
- Surrogate boundary: `integrated_latent.h5mu` is a MuData-inspired HDF5 layout written with `h5py`, not a full `muon` serialization.

## Package scaffold

- `metadata.yaml` defines the experiment contract and required deliverables.
- `examples/toy_input.json` stores raw synthetic multimodal counts and expected invariants only.
- `scripts/run_multimodal_integration.py` runs the deterministic starter on raw toy inputs.
- `scripts/validate_outputs.py` checks the deliverables and the MuData-inspired HDF5 structure.
- `tests/test_multimodal_integration.py` asserts numeric and structural invariants of the starter path.

## Toy validation

- Run `python3 scripts/run_exercise.py --outdir scratch/multimodal-integration_toy_run`.
- The example finishes in a few seconds on eight synthetic cells across two batches.
- Treat the outputs as starter validation, not biological evidence.

## Real-run expectations

- Replace the toy counts with project or public inputs that preserve the same deliverable names.
- Keep model choice, modality weights, feature definitions, and label-transfer references explicit in the final report.
- Stop and document blockers if identifiers, batch metadata, or reference labels are missing.

## Failure handling

- Fail early when cell identifiers are duplicated, modality dimensions do not align, or reference/query splits are incomplete.
- Preserve uncertainty in the report instead of overstating biological certainty.
- Do not treat the starter as a substitute for full totalVI or MultiVI training on real multimodal datasets.
