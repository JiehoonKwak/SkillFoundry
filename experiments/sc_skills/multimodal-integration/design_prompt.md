Design or refine an experiment-only portable skill for the task: Multimodal Integration.

Task goal:
Integrate RNA, protein, and chromatin modalities with explicit model choice, latent-space QC, and label-transfer outputs.

Constraints:
- Only create or modify files under `experiments/sc_skills/multimodal-integration`.
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
- totalVI or MultiVI
- latent-space QC
- label transfer

Family-specific starter guidance:
- Do not require true totalVI or MultiVI training; implement a lightweight surrogate starter that still computes something real.
- Starter operations should normalize RNA and protein separately, concatenate or use a tiny CCA/PCA/SVD latent, build a kNN graph, run label transfer from reference half to query half, and emit batch-mixing plus label-consistency QC.
- Keep the toy problem very small, such as 8 cells, 4 RNA genes, 2 proteins, and 2 batches.
- The latent representation and QC must be computed from raw toy inputs, not stored directly in a pseudo `.h5mu`.

Required deliverables:
- integrated_latent.h5mu
- modality_qc.tsv
- integration_report.md

Refresh the experiment package to keep its docs, scripts, examples, and tests aligned.
Add or refine a tiny runnable example so `scripts/run_exercise.py` performs a meaningful method-shaped starter path instead of only echoing precomputed outputs.
Update `refs.md` if you need additional official docs, papers, notebooks, or GitHub references to justify the starter design.
