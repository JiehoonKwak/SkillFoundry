---
name: database-query-portable-skill
description: Use this experiment skill to answer small gene, protein, variant, and disease lookup tasks with deterministic identifier normalization, mock public-API fan-out, and cross-source reconciliation that stays portable across Codex, Claude Code, and similar shell-capable agents.
---
## Purpose

Answer small bioinformatics database-query tasks by normalizing raw query strings, searching bundled alias tables plus mock public-source payloads, and reconciling the hits into stable identifiers and reviewable notes.

## When to use

- You need a portable starter for gene, protein, variant, or disease lookup logic before wiring in live APIs.
- You want a deterministic local exercise path that mirrors Ensembl, MyGene.info, UniProt, ClinVar, MedGen, and OMIM style lookups without network access.
- You need machine-readable QC for normalization, fan-out, and conflict resolution.

## When not to use

- You need authoritative live annotations, current review status, or full database coverage.
- You need full HGVS parsing, transcript-aware variant projection, or clinical-grade interpretation.
- You need to query private mirrors, credentials, or large local parquet stores from the original Spatial Agent setup.

## Inputs

- `examples/toy_input.json` or a compatible JSON file with:
- raw query strings,
- local concept and alias tables,
- mock source payloads,
- source-priority weights,
- expected invariants for the toy exercise.

## Outputs

- `query_results.tsv`
- `source_provenance.md`
- `resolution_notes.md`
- `normalization_qc.tsv`
- `fanout_qc.tsv`
- `reconciliation_qc.tsv`
- `run_summary.json`

## Requirements

- `python3`
- No network access
- No private Spatial Agent helpers

## Procedure

1. Normalize each raw query into reusable lookup keys.
2. Infer whether the query behaves like a gene, protein, variant, or disease request.
3. Fan out over bundled alias tables and mock Ensembl, MyGene.info, UniProt, ClinVar, MedGen, and OMIM payloads.
4. Reconcile competing concepts with deterministic source weights and preferred-ID rules.
5. Write the declared deliverables plus machine-readable QC artifacts.

## Validation

- `python3 scripts/run_database_query.py --input examples/toy_input.json --outdir scratch/database-query_toy_run`
- `python3 scripts/validate_outputs.py --input examples/toy_input.json --outdir scratch/database-query_toy_run`
- `python3 -m unittest tests/test_contract.py tests/test_database_query.py`

## Starter scope

- Truly computed: query normalization, entity-type scoring, source fan-out over bundled mock payloads, and score-based reconciliation into normalized identifiers.
- Approximated: source ranking heuristics, simplified payload schemas, and local crosswalks that stand in for richer live API joins.
- Surrogate only: live HTTP requests, full HGVS validation, transcript selection, pagination, retry logic, and clinical evidence review.

## Failure modes and fixes

- If no concept is found, add aliases or source records instead of hard-coding deliverables.
- If a query resolves to the wrong entity type, tune the query-type keywords or source weights in the input JSON.
- If two concepts stay close after scoring, keep the ambiguity explicit in `resolution_notes.md` and inspect `reconciliation_qc.tsv`.

## Safety and limits

- Treat the starter outputs as workflow validation, not biological or clinical evidence.
- Do not use the toy run to justify patient-care, wet-lab, or regulatory decisions.
- Replace the mock payloads with pinned live-source snapshots before promoting this experiment toward production use.

## Examples

- Run the bundled toy exercise:
- `python3 scripts/run_exercise.py --outdir scratch/toy_run`
- Swap in a new JSON payload with the same top-level structure to test different aliases or scoring rules.

## Provenance

- Adapted from `source_reference.md` for the same task family.
- Private Spatial Agent path assumptions and hidden database wrappers were removed.
- Public method anchors are tracked in `refs.md`.

## Related skills

- None inside `experiments/sc_skills/`; this package is a task-specific micro-starter for database lookup logic.
