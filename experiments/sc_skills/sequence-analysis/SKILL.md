---
name: sequence-analysis-portable-skill
description: Use this experiment skill to identify short DNA sequences, annotate overlapping features, inspect tracked variants, and rank primer-oriented follow-up candidates in one portable workflow for Codex, Claude Code, and similar shell-capable agents.
---
## Purpose

Handle compact sequence-identification and primer-design tasks with a deterministic local starter that stays runnable without BLAST databases, Primer3 binaries, or network access.

## When to use

- You need a small portable starter for sequence lookup, feature overlap, variant inspection, and primer follow-up on bundled toy inputs.
- You want machine-readable QC for alignment ranking, reverse-complement handling, tracked variant calls, and primer-pair scoring.
- You need an experiment-only scaffold that can later be swapped to BLAST+, Ensembl, or Primer3-backed runs while keeping the same deliverable names.

## When not to use

- You need genome-scale alignments, transcript-aware annotations, or current public database answers.
- You need full Primer3 thermodynamic modeling, mispriming search against a genome, or clinical-grade variant interpretation.
- You need wet-lab or patient-care conclusions from the toy outputs.

## Inputs

- `examples/toy_input.json` or a compatible JSON file with:
- short query sequences,
- bundled reference slices and interval annotations,
- tracked variant positions and expected alleles,
- lookup thresholds and primer constraints,
- expected invariants for the toy exercise.

## Outputs

- `sequence_summary.tsv`
- `feature_annotations.tsv`
- `primer_candidates.tsv`
- `lookup_qc.tsv`
- `variant_qc.tsv`
- `primer_qc.tsv`
- `run_summary.json`

## Requirements

- `python3`
- No network access
- No private Spatial Agent helpers

## Procedure

1. Normalize each query sequence and score ungapped forward and reverse-complement alignments against the bundled references.
2. Keep the top-scoring reference hit, record lookup QC, and annotate overlapping features plus tracked variant sites.
3. Call the tracked variant state directly from the aligned toy sequence.
4. Scan the aligned query for primer pairs that satisfy GC, Wallace Tm, GC-clamp, spacing, and product-size rules, then rank them with simple penalties.
5. Validate the declared deliverables and starter invariants.

## Validation

- `python3 scripts/run_sequence_analysis.py --input examples/toy_input.json --outdir scratch/sequence-analysis_toy_run`
- `python3 scripts/validate_outputs.py --input examples/toy_input.json --outdir scratch/sequence-analysis_toy_run`
- `python3 -m unittest tests/test_contract.py tests/test_sequence_analysis.py`

## Starter scope

- Truly computed: ungapped sequence lookup, reverse-complement handling, interval overlap, tracked variant-state inspection, and ranked primer scans on the synthetic slices.
- Approximated: BLAST-style search is replaced by windowed overlap and longest-seed scoring; Primer3 is replaced by Wallace Tm, GC%, GC-clamp, homopolymer, and simple self-complement heuristics.
- Surrogate only: live BLAST databases, Primer3 thermodynamic models, transcript-aware feature retrieval, off-target genome searches, and real clinical variant interpretation.

## Failure modes and fixes

- If the wrong reference wins, adjust the bundled sequences or `lookup_policy` thresholds instead of hard-coding the summary table.
- If no primer pairs pass, widen `product_size_range`, relax the GC/Tm bounds slightly, or provide a longer query slice around the target site.
- If a tracked variant is outside the aligned region, extend the query sequence or move the tracked position in the toy reference so the exercise remains explicit.

## Safety and limits

- Treat the starter outputs as workflow validation on synthetic sequences, not biological evidence.
- Do not use the toy run for experimental design claims beyond smoke-testing the portable contract.
- Replace the toy slices with pinned public references and real tool calls before promoting this experiment toward production use.

## Examples

- Run the bundled toy exercise:
- `python3 scripts/run_exercise.py --outdir scratch/toy_run`
- Swap in a new JSON payload with the same top-level structure to test a different synthetic hotspot or primer constraint set.

## Provenance

- Adapted from `source_reference.md` as decomposition guidance only.
- Private Spatial Agent dependencies were removed from the runnable path.
- Public method anchors for BLAST, Primer3, and sequence-feature services are tracked in `refs.md`.

## Related skills

- None inside `experiments/sc_skills/`; this package is the task-specific micro-starter for portable sequence lookup and primer follow-up.
