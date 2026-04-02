Design or refine an experiment-only portable skill for the task: Spatial Panel Design.

Task goal:
Design a balanced marker panel for spatial transcriptomics or targeted assays using explicit evidence and redundancy control.

Constraints:
- Only create or modify files under `experiments/sc_skills/panel-design`.
- Do not edit `registry/`, `skills/`, `site/`, `README.md`, `experiments.md`, or planning files.
- Keep the result compatible with Codex, Claude Code, and similar shell-capable agents.
- Use `source_reference.md` as decomposition guidance only; remove private Spatial Agent tool dependencies.
- Keep `metadata.yaml`, `examples/`, `scripts/`, `tests/`, and `assets/` aligned with the experiment contract.
- Upgrade the package from a static toy contract toward a small runnable starter that still runs quickly on tiny synthetic or pseudo inputs.
- Prefer method-shaped computations, heuristics, rankings, graph operations, or QC logic over directly dumping fixed files whenever that is feasible on toy data.
- Preserve the declared deliverable names and schemas so repository-level validation still passes.
- If a real public dependency is too heavy for the tiny starter, replace it with a documented lightweight approximation and state the boundary clearly in `SKILL.md` and `assets/README.md`.

Hard requirements for this refine pass:
- Do not use `examples/toy_input.json` to store final deliverables or narrative report text. The toy input may only hold raw synthetic inputs, priors, reference tables, coordinates, or expected invariants.
- Replace the shared contract-only runner with a task-specific runner that computes outputs from the raw toy inputs using deterministic local operations and no network access.
- Implement 2 to 4 method-shaped computations that finish in under 5 seconds on tiny synthetic inputs and emit machine-readable intermediate QC, not only final deliverables.
- Tests must assert numeric or structural invariants of the computation, not only file existence.
- Add a short `Starter scope` note to `SKILL.md` stating what is truly computed, what is approximated, and what remains a surrogate for the full public method.

Preferred methods:
- marker aggregation
- redundancy pruning
- coverage balancing

Family-specific starter guidance:
- Do not hide multiple unrelated tasks behind one abstract template; each task should get a task-specific micro-starter.
- `database-query`: use raw toy inputs like alias tables, mock API JSON payloads, and query lists; compute identifier normalization, source merge, and conflict resolution.
- `sequence-analysis`: use raw toy inputs like short FASTA records, feature intervals, and variant positions; compute exact or approximate lookup, reverse complement, and a simple primer scan with GC/Tm/product-size rules.
- `panel-design`: use raw toy inputs like marker score tables and cell-type coverage matrices; compute marker ranking, redundancy pruning, and coverage balancing.

Required deliverables:
- panel_candidates.tsv
- panel_rationale.md
- platform_notes.md

Refresh the experiment package to keep its docs, scripts, examples, and tests aligned.
Add or refine a tiny runnable example so `scripts/run_exercise.py` performs a meaningful method-shaped starter path instead of only echoing precomputed outputs.
Update `refs.md` if you need additional official docs, papers, notebooks, or GitHub references to justify the starter design.
