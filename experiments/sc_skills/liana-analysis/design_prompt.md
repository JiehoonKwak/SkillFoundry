Design or refine an experiment-only portable skill for the task: LIANA Communication Analysis.

Task goal:
Run LIANA or LIANA+ with explicit method selection, harmonized metadata, and interpretable interaction prioritization.

Constraints:
- Only create or modify files under `experiments/sc_skills/liana-analysis`.
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
- LIANA ranking
- resource harmonization
- condition-aware summaries

Family-specific starter guidance:
- Use raw toy inputs that look like cell or group expression matrices, group labels, a ligand-receptor catalog, and an optional adjacency matrix.
- Starter operations should include group means, ligand-receptor scoring with `min(ligand_expr, receptor_expr)` or a geometric mean, a small shuffle-based null ranking, and optional adjacency support.
- Keep the toy problem very small, such as 8 cells, 3 cell groups, 4 ligand-receptor pairs, and one known positive pair.
- For CellPhoneDB or LIANA style tasks, do not pretend to run the full heavy stack; at minimum compute expression filtering, ranking, and a tiny null model.

Required deliverables:
- liana_rankings.tsv
- resource_overlap.tsv
- liana_report.md

Refresh the experiment package to keep its docs, scripts, examples, and tests aligned.
Add or refine a tiny runnable example so `scripts/run_exercise.py` performs a meaningful method-shaped starter path instead of only echoing precomputed outputs.
Update `refs.md` if you need additional official docs, papers, notebooks, or GitHub references to justify the starter design.
