# Coverage Report

## Snapshot

- Top-level taxonomy domains: `27`
- Canonical resources in `resources_dedup.jsonl`: `375`
- Skills in `skills.jsonl`: `267`
- `sandbox_verified` skills: `118`
- `slurm_verified` skills: `4`
- Covered domains in `site/tree.json`: `27`
- Covered leaves: `254`
- Frontier-only leaves: `0`
- Pure TODO leaves: `0`

## Phase 63

Phase 63 added a concrete local protein language model starter for embedding extraction and sequence-to-function triage, bridging the previously unverified proteomics frontier placeholders with a reusable executable contract.

Key findings:

- The registry expanded to `375` canonical resources and `267` skills.
- `protein-language-model-function-analysis-starter` validates protein FASTA input, emits a tabular embedding matrix, computes nearest-neighbor structure, and runs centroid-based label transfer on toy examples.
- The new skill is `sandbox_verified`, increasing the repository-wide `sandbox_verified` total to `118`.
- Coverage totals stayed unchanged at `254` covered leaves, `0` frontier-only leaves, and `0` pure TODO leaves because this pass deepened an already covered proteomics area.

Validation for this phase:

- `python3 skills/proteomics/protein-language-model-function-analysis-starter/scripts/run_protein_language_model_function_analysis.py --input skills/proteomics/protein-language-model-function-analysis-starter/examples/toy_sequences.fasta --labels skills/proteomics/protein-language-model-function-analysis-starter/examples/toy_labels.tsv --config skills/proteomics/protein-language-model-function-analysis-starter/examples/analysis_config.json --embeddings-out scratch/protein-lm/toy_embeddings.tsv --summary-out scratch/protein-lm/toy_summary.json`
- `python3 -m unittest discover -s skills/proteomics/protein-language-model-function-analysis-starter/tests -p 'test_*.py'`
- `python3 -m unittest tests.smoke.test_protein_language_model_function_analysis_starter`
- `python3 skills/proteomics/protein-language-model-function-analysis-starter/scripts/run_frontier_starter.py --out scratch/protein-lm/frontier_summary.json`
- `make smoke-protein-language-models`
- `python3 scripts/validate_repository.py`
- `python3 scripts/build_site.py`

## Phase 62

Phase 62 upgraded the existing `gwas-starter` from a frontier placeholder into a runnable GWAS summary-statistics QC skill with deterministic example data and downstream interpretation guidance.

Key findings:

- The registry expanded to `371` canonical resources while the skill count stayed at `266`.
- `gwas-starter` now exposes a stable local contract: read a GWAS sumstats table, normalize common headers, flag duplicate or malformed variants, write a flagged TSV, and summarize next steps for PLINK 2.0, LDSC, FUMA, GWASLab, and GWAS Catalog SSF follow-up.
- The promoted skill moved from `implemented` to `sandbox_verified`, increasing the repository-wide `sandbox_verified` total to `117`.
- The site-level coverage counts stayed unchanged at `254` covered leaves, `0` frontier-only leaves, and `0` pure TODO leaves because this pass refined an already covered leaf instead of adding a new one.

Validation for this phase:

- `python3 skills/genomics/gwas-starter/scripts/run_gwas_summary_qc.py --input skills/genomics/gwas-starter/examples/toy_sumstats.tsv --config skills/genomics/gwas-starter/examples/qc_config.json --out-tsv skills/genomics/gwas-starter/assets/toy_gwas_qc.tsv --summary-out skills/genomics/gwas-starter/assets/toy_gwas_qc_summary.json`
- `python3 -m unittest discover -s skills/genomics/gwas-starter/tests -p 'test_*.py'`
- `python3 scripts/validate_repository.py`
- `python3 skills/genomics/gwas-starter/scripts/run_frontier_starter.py --out scratch/gwas/frontier_wrapper_summary.json`
- `python3 scripts/build_site.py`
- `python3 -m unittest tests.integration.test_build_site`

## Phase 61

Phase 61 was a conservative `refresh` pass that rebuilt the generated site files, refreshed the shared review outputs, and wrote a new dated per-run report bundle without changing registry truth or verification labels.

This pass did not expand the taxonomy or registries. Its job was to keep the refresh-owned artifacts internally consistent after the current framework and repository changes already present in the checkout.

Key findings:

- The repository still sits at `366` canonical resources, `266` skills, `254` covered leaves, `0` frontier-only leaves, and `0` pure TODO leaves.
- `scripts/close_frontier_leaves.py` remained a no-op, which is still the correct behavior for the fully skill-backed taxonomy.
- `scripts/audit_skill_suite.py` stayed at `0` hard failures, `266/266` smoke-covered skills, and `253/266` skills with local tests.
- `scripts/run_skill_smoke_matrix.py --dry-run` still resolved `266/266` smoke targets with no missing mappings, and this phase kept that result labeled as mapping-only verification.
- The targeted six-case benchmark suite covering GWAS, proteomics, privacy-preserving analysis, pseudobulk analysis, and interface analysis kept maintained wrappers at `1.0` average deliverable rate versus `0.611` for the baseline while both paths still succeeded at command execution.
- `python3 -m unittest tests.integration.test_build_site tests.integration.test_skill_review_scripts tests.integration.test_sciskill_framework` passed `52` tests, keeping the framework-facing refresh surface conservative and green.
- `scripts/build_site.py` kept the generated tree totals unchanged while refreshing the generated site outputs and current framework-run index.

Refined carry-forward queue for the requested priorities:

- `proteomics :: differential-proteomics`
- `proteomics :: ms-proteomics-preprocessing`
- `proteomics :: sequence-to-function-modeling`
- `proteomics :: protein-embeddings`
- `proteomics :: interface-analysis`
- `transcriptomics :: spatial-transcriptomics`
- `transcriptomics :: pseudobulk-analysis`
- `genomics :: gwas`
- `clinical_biomedical_data_science :: privacy-preserving-analysis`
- `systems_biology_and_network_science :: regulatory-network-inference`

Validation for this phase:

- `python3 scripts/validate_repository.py`
- `python3 scripts/build_site.py`
- `python3 scripts/close_frontier_leaves.py`
- `python3 scripts/audit_skill_suite.py --json-out scratch/refresh_audit.json --markdown-out scratch/refresh_audit.md`
- `python3 scripts/run_skill_smoke_matrix.py --dry-run --json-out scratch/refresh_smoke_matrix.json --markdown-out scratch/refresh_smoke_matrix.md`
- `python3 scripts/benchmark_skill_advantage.py --case gwas-starter-summary --case ms-proteomics-preprocessing-starter-summary --case protein-embeddings-starter-summary --case privacy-preserving-analysis-starter-summary --case pseudobulk-analysis-starter-summary --case interface-analysis-starter-summary --json-out scratch/refresh_benchmark.json --markdown-out scratch/refresh_benchmark.md`
- `python3 -m unittest tests.integration.test_build_site tests.integration.test_skill_review_scripts tests.integration.test_sciskill_framework`

## Phase 60

Phase 60 was a conservative `refresh` pass that refreshed the shared reports and generated outputs after the latest framework and documentation changes already present in the checkout.

This pass did not expand the taxonomy or registries. Its job was to keep the refresh-owned artifacts internally consistent while recording exactly which checks were mapping-only versus actually executed.

Key findings:

- The repository still sits at `366` canonical resources, `266` skills, `254` covered leaves, `0` frontier-only leaves, and `0` pure TODO leaves.
- `scripts/close_frontier_leaves.py` remained a no-op, which is still the correct behavior for the fully skill-backed taxonomy.
- `scripts/audit_skill_suite.py` stayed at `0` hard failures, `266/266` smoke-covered skills, `253/266` skills with local tests, and `13` warning-only gaps.
- `scripts/run_skill_smoke_matrix.py --dry-run` still resolved `266/266` smoke targets with no missing mappings, but this remained a mapping refresh rather than a full suite execution.
- The focus-leaf execution lane ran the `12` surfaced starter skills for differential proteomics, MS proteomics preprocessing, sequence-to-function modeling, climate reanalysis access, physics-informed forecasting, raster/vector ingestion, scientific map/dashboard generation, spatial interpolation, GWAS, symbolic discovery of physical laws, interface analysis, and protein embeddings.
- The targeted benchmark suite covered `12` representative cases for GWAS, proteomics, climate reanalysis, single-cell planning, and multiome integration; maintained skills kept equal success with a higher mean deliverable rate (`1.0` versus `0.583`).
- `scripts/build_site.py` kept the generated tree totals unchanged while surfacing `16` framework runs in `site/framework_runs.json`.

Refined carry-forward queue for the requested priorities:

- `proteomics :: differential-proteomics`
- `proteomics :: ms-proteomics-preprocessing`
- `proteomics :: sequence-to-function-modeling`
- `proteomics :: protein-embeddings`
- `proteomics :: interface-analysis`
- `transcriptomics :: spatial-transcriptomics`
- `transcriptomics :: pseudobulk-analysis`
- `transcriptomics :: rna-velocity`
- `genomics :: gwas`
- `genomics :: fine-mapping`
- `clinical_biomedical_data_science :: privacy-preserving-analysis`
- `systems_biology_and_network_science :: regulatory-network-inference`

Validation for this phase:

- `python3 scripts/close_frontier_leaves.py`
- `python3 scripts/validate_repository.py`
- `python3 scripts/audit_skill_suite.py --json-out reports/refresh-runs/20260321-050747-cycle/audit_skill_suite.json --markdown-out reports/refresh-runs/20260321-050747-cycle/audit_skill_suite.md`
- `python3 scripts/run_skill_smoke_matrix.py --dry-run --json-out reports/refresh-runs/20260321-050747-cycle/run_skill_smoke_matrix.json --markdown-out reports/refresh-runs/20260321-050747-cycle/run_skill_smoke_matrix.md`
- Focus smoke commands from the `12` surfaced starter skills, logged in `reports/refresh-runs/20260321-050747-cycle/focus_smoke.log`
- `python3 scripts/benchmark_skill_advantage.py --case gwas-starter-summary --case gwas-starter-checklist --case interface-analysis-starter-summary --case interface-analysis-starter-checklist --case ms-proteomics-preprocessing-starter-summary --case ms-proteomics-preprocessing-starter-checklist --case climate-reanalysis-access-starter-summary --case climate-reanalysis-access-starter-checklist --case langgraph-planning-execution-agent-single-cell-report --case langgraph-planning-execution-agent-literature-single-cell-report --case multiome-integration-starter-summary --case multiome-integration-starter-augmented --json-out reports/refresh-runs/20260321-050747-cycle/benchmark_skill_advantage.json --markdown-out reports/refresh-runs/20260321-050747-cycle/benchmark_skill_advantage.md`
- `python3 scripts/build_site.py`

## Phase 59

Phase 59 was a conservative `refresh` pass that reran the shared repository review surface, refreshed the per-run artifacts, and rebuilt the generated site outputs after the current framework changes.

This pass did not change the taxonomy, resource registry, or skill registry. Its job was to confirm that the refresh-owned shared files still agree with the registries and that the framework remains consistent after rebuilding the generated outputs.

Key findings:

- The repository remained at `366` canonical resources, `266` skills, `254` covered leaves, `0` frontier-only leaves, and `0` pure TODO leaves.
- `scripts/close_frontier_leaves.py` remained a no-op, which is still the correct behavior for the fully skill-backed taxonomy.
- `scripts/audit_skill_suite.py` stayed at `0` hard failures, `266/266` smoke-covered skills, `253/266` skills with local tests, and `13` warning-only gaps.
- `scripts/run_skill_smoke_matrix.py --dry-run` still resolved `266/266` smoke targets with no missing mappings.
- The focused `fastqc-multiqc-read-qc` benchmark still showed equal command success for skill and baseline while preserving the maintained skill's higher deliverable rate (`1.0` versus `0.75`).
- `python3 -m unittest tests.integration.test_sciskill_framework tests.integration.test_skill_review_scripts -v` passed all `22` selected integration tests, covering the campaign-runner and hierarchical-evaluation framework paths.
- `scripts/build_site.py` kept the generated tree totals unchanged while surfacing `12` framework runs in `site/framework_runs.json`, all still counted conservatively under `attention`.

Refined carry-forward queue for the requested priorities:

- `proteomics :: differential-proteomics`
- `proteomics :: ms-proteomics-preprocessing`
- `proteomics :: sequence-to-function-modeling`
- `proteomics :: protein-embeddings`
- `proteomics :: protein-identification-quantification`
- `transcriptomics :: spatial-transcriptomics`
- `transcriptomics :: pseudobulk-analysis`
- `transcriptomics :: rna-velocity`
- `genomics :: gwas`
- `genomics :: fine-mapping`
- `clinical_biomedical_data_science :: privacy-preserving-analysis`
- `systems_biology_and_network_science :: regulatory-network-inference`

Validation for this phase:

- `python3 scripts/validate_repository.py`
- `python3 scripts/close_frontier_leaves.py`
- `python3 scripts/audit_skill_suite.py --json-out reports/refresh-20260321-cycle/audit_skill_suite.json --markdown-out reports/refresh-20260321-cycle/audit_skill_suite.md`
- `python3 scripts/run_skill_smoke_matrix.py --dry-run --json-out reports/refresh-20260321-cycle/smoke_matrix_dry_run.json --markdown-out reports/refresh-20260321-cycle/smoke_matrix_dry_run.md`
- `python3 scripts/benchmark_skill_advantage.py --case fastqc-multiqc-read-qc --json-out reports/refresh-20260321-cycle/benchmark_fastqc.json --markdown-out reports/refresh-20260321-cycle/benchmark_fastqc.md`
- `python3 -m unittest tests.integration.test_sciskill_framework tests.integration.test_skill_review_scripts -v`
- `python3 scripts/build_site.py`

## Phase 58

Phase 58 was a conservative `tree_check` pass that refreshed the shared audit surface, rebuilt the site outputs, and re-ranked the next verification queue around the requested protein, spatial omics, sequence analysis, GWAS, and privacy or regulatory priorities.

This pass did not change the taxonomy, resource registry, or skill registry. Its job was to confirm that the zero-frontier tree remains internally consistent and that the next backlog is still the unverified implementation queue rather than breadth expansion.

Key findings:

- The repository remained at `366` canonical resources, `266` skills, `254` covered leaves, `0` frontier-only leaves, and `0` pure TODO leaves.
- `scripts/close_frontier_leaves.py` remained a no-op, which is still the correct behavior for a fully skill-backed taxonomy.
- `scripts/audit_skill_suite.py` stayed at `0` hard failures, `266/266` smoke-covered skills, `253/266` skills with local tests, and `13` warning-only gaps.
- `scripts/run_skill_smoke_matrix.py --dry-run` still resolved `266/266` smoke targets with no missing mappings.
- The focused `fastqc-multiqc-read-qc` benchmark still showed equal command success for skill and baseline while preserving the maintained skill's higher deliverable rate (`1.0` versus `0.75`).
- A refreshed `scripts/sciskill_framework.py status` snapshot confirmed that `single-cell-rna-seq-preprocessing` and `single-cell-integration-batch-correction` already have verified skills, so the more valuable transcriptomics backlog remains `spatial-transcriptomics`, `pseudobulk-analysis`, and `rna-velocity`.
- The broad `legal` and `policy` terms still collapse onto `privacy-preserving-analysis`, so the recommended follow-up is to keep future focus terms explicit and conservative instead of changing framework behavior during tree-check.

Refined carry-forward queue for the requested priorities:

- `proteomics :: differential-proteomics`
- `proteomics :: ms-proteomics-preprocessing`
- `proteomics :: sequence-to-function-modeling`
- `proteomics :: protein-embeddings`
- `proteomics :: protein-identification-quantification`
- `transcriptomics :: spatial-transcriptomics`
- `transcriptomics :: pseudobulk-analysis`
- `transcriptomics :: rna-velocity`
- `genomics :: gwas`
- `genomics :: fine-mapping`
- `clinical_biomedical_data_science :: privacy-preserving-analysis`
- `systems_biology_and_network_science :: regulatory-network-inference`

Validation for this phase:

- `python3 scripts/validate_repository.py`
- `python3 scripts/close_frontier_leaves.py`
- `python3 scripts/audit_skill_suite.py --json-out reports/tree-check-20260321-cycle/audit_skill_suite.json --markdown-out reports/tree-check-20260321-cycle/audit_skill_suite.md`
- `python3 scripts/run_skill_smoke_matrix.py --dry-run --json-out reports/tree-check-20260321-cycle/smoke_matrix_dry_run.json --markdown-out reports/tree-check-20260321-cycle/smoke_matrix_dry_run.md`
- `python3 scripts/benchmark_skill_advantage.py --case fastqc-multiqc-read-qc --json-out reports/tree-check-20260321-cycle/benchmark_fastqc.json --markdown-out reports/tree-check-20260321-cycle/benchmark_fastqc.md`
- `python3 scripts/sciskill_framework.py --json status --focus-limit 24 --focus-term proteomics --focus-term protein --focus-term sequence --focus-term gwas --focus-term "single-cell" --focus-term "single cell" --focus-term pseudobulk --focus-term "rna velocity" --focus-term "spatial transcriptomics" --focus-term privacy --focus-term regulatory --focus-term legal --focus-term policy > reports/tree-check-20260321-cycle/framework_status.json`
- `python3 scripts/build_site.py`

## Phase 57

Phase 57 was a conservative `refresh` pass that reran the shared repository review surface, refreshed the per-run artifacts, and revalidated the framework-facing integration path after the recent framework changes.

This pass did not change the taxonomy, resource registry, or skill registry. Its job was to confirm that the refresh-owned shared files still agree with the registries and that the framework remains consistent after rebuilding the generated site outputs.

Key findings:

- The repository remained at `366` canonical resources, `266` skills, `254` covered leaves, `0` frontier-only leaves, and `0` pure TODO leaves.
- `scripts/close_frontier_leaves.py` remained a no-op, which is still the correct behavior for the fully skill-backed taxonomy.
- `scripts/audit_skill_suite.py` stayed at `0` hard failures, `266/266` smoke-covered skills, `253/266` skills with local tests, and `13` warning-only gaps.
- `scripts/run_skill_smoke_matrix.py --dry-run` still resolved `266/266` smoke targets with no missing mappings.
- The focused `fastqc-multiqc-read-qc` benchmark still showed equal command success for skill and baseline while preserving the maintained skill's higher deliverable rate (`1.0` versus `0.75`).
- `python3 -m unittest tests.integration.test_sciskill_framework tests.integration.test_skill_review_scripts -v` passed all `17` integration tests, so the refresh surface remains framework-safe after the recent orchestration changes.
- `scripts/build_site.py` kept the generated site totals unchanged while surfacing `8` framework runs in `site/framework_runs.json`.

Refined carry-forward queue for the requested priorities:

- `proteomics :: differential-proteomics`
- `proteomics :: ms-proteomics-preprocessing`
- `proteomics :: sequence-to-function-modeling`
- `proteomics :: protein-embeddings`
- `proteomics :: protein-identification-quantification`
- `transcriptomics :: spatial-transcriptomics`
- `transcriptomics :: pseudobulk-analysis`
- `transcriptomics :: rna-velocity`
- `genomics :: gwas`
- `genomics :: fine-mapping`
- `clinical_biomedical_data_science :: privacy-preserving-analysis`
- `systems_biology_and_network_science :: regulatory-network-inference`

Validation for this phase:

- `python3 scripts/validate_repository.py`
- `python3 scripts/close_frontier_leaves.py`
- `python3 scripts/audit_skill_suite.py --json-out reports/refresh-20260321-cycle/audit_skill_suite.json --markdown-out reports/refresh-20260321-cycle/audit_skill_suite.md`
- `python3 scripts/run_skill_smoke_matrix.py --dry-run --json-out reports/refresh-20260321-cycle/smoke_matrix_dry_run.json --markdown-out reports/refresh-20260321-cycle/smoke_matrix_dry_run.md`
- `python3 scripts/benchmark_skill_advantage.py --case fastqc-multiqc-read-qc --json-out reports/refresh-20260321-cycle/benchmark_fastqc.json --markdown-out reports/refresh-20260321-cycle/benchmark_fastqc.md`
- `python3 -m unittest tests.integration.test_sciskill_framework tests.integration.test_skill_review_scripts -v`
- `python3 scripts/build_site.py`

## Phase 56

Phase 56 was a conservative `tree_check` pass that refreshed the audit surface after the recent framework work and narrowed the next implementation queue to the requested scientific priorities.

This pass did not change the taxonomy, resource registry, or skill registry. Its job was to confirm that the tree remains fully covered and to distinguish already-verified single-cell leaves from the still-unverified protein, transcriptomics, GWAS, and privacy or regulatory starters.

Key findings:

- The repository remained at `366` canonical resources, `266` skills, `254` covered leaves, `0` frontier-only leaves, and `0` pure TODO leaves.
- `scripts/close_frontier_leaves.py` remained a no-op, which is still the correct behavior for a fully skill-backed taxonomy.
- `scripts/audit_skill_suite.py` stayed at `0` hard failures, `266/266` smoke-covered skills, `253/266` skills with local tests, and `13` warning-only gaps.
- `scripts/run_skill_smoke_matrix.py --dry-run` still resolved `266/266` smoke targets with no missing mappings.
- The focused `fastqc-multiqc-read-qc` benchmark still showed equal command success for skill and baseline while preserving the maintained skill's stronger deliverable structure.
- A narrower `scripts/sciskill_framework.py status` query showed that `single-cell-rna-seq-preprocessing` and `single-cell-integration-batch-correction` already have verified skills, so the more valuable next-pass transcriptomics targets are `spatial-transcriptomics`, `pseudobulk-analysis`, and `rna-velocity`.

Refined carry-forward queue for the requested priorities:

- `proteomics :: differential-proteomics`
- `proteomics :: ms-proteomics-preprocessing`
- `proteomics :: sequence-to-function-modeling`
- `proteomics :: protein-embeddings`
- `proteomics :: protein-identification-quantification`
- `transcriptomics :: spatial-transcriptomics`
- `transcriptomics :: pseudobulk-analysis`
- `transcriptomics :: rna-velocity`
- `genomics :: gwas`
- `genomics :: fine-mapping`
- `clinical_biomedical_data_science :: privacy-preserving-analysis`
- `systems_biology_and_network_science :: regulatory-network-inference`

Validation for this phase:

- `python3 scripts/validate_repository.py`
- `python3 scripts/close_frontier_leaves.py`
- `python3 scripts/audit_skill_suite.py --json-out reports/tree-check-20260321-cycle/audit_skill_suite.json --markdown-out reports/tree-check-20260321-cycle/audit_skill_suite.md`
- `python3 scripts/run_skill_smoke_matrix.py --dry-run --json-out reports/tree-check-20260321-cycle/smoke_matrix_dry_run.json --markdown-out reports/tree-check-20260321-cycle/smoke_matrix_dry_run.md`
- `python3 scripts/benchmark_skill_advantage.py --case fastqc-multiqc-read-qc --json-out reports/tree-check-20260321-cycle/benchmark_fastqc.json --markdown-out reports/tree-check-20260321-cycle/benchmark_fastqc.md`
- `python3 scripts/sciskill_framework.py --json status --focus-limit 20 --focus-term proteomics --focus-term protein --focus-term sequence --focus-term gwas --focus-term single-cell --focus-term "single cell" --focus-term pseudobulk --focus-term rna-velocity --focus-term spatial-transcriptomics --focus-term privacy --focus-term regulatory`
- `python3 scripts/build_site.py`

## Phase 55

Phase 55 was a conservative `refresh` pass that rebuilt the generated site artifacts and refreshed the report surfaces after the earlier framework changes.

This pass did not change the taxonomy, resource registry, or skill registry. Its job was to confirm that the repository-generated outputs still agree with the registries and to write a fresh per-run summary.

Key findings:

- The repository remained at `366` canonical resources, `266` skills, `254` covered leaves, `0` frontier-only leaves, and `0` pure TODO leaves.
- `scripts/close_frontier_leaves.py` stayed a no-op, which is the correct behavior now that every taxonomy leaf is already skill-backed.
- `scripts/audit_skill_suite.py` stayed at `0` hard failures, `266/266` smoke-covered skills, `253/266` skills with local tests, and `13` warning-only gaps.
- `scripts/run_skill_smoke_matrix.py --dry-run` still resolved `266/266` smoke targets with no missing mappings.
- The focused `fastqc-multiqc-read-qc` benchmark still showed equal command success for skill and baseline while preserving the maintained skill's stronger deliverable rate (`1.0` versus `0.75`).
- `scripts/build_site.py` refreshed `site/framework_runs.json` to include `6` manifest-backed framework runs plus the current `latest_status` snapshot from `scratch/framework/latest_status.json`.
- All surfaced framework runs currently remain in `attention` state, which matches the stored manifests and keeps the website view conservative rather than silently downgrading outstanding parser or novelty-review issues.

Refined carry-forward queue for the requested priorities:

- `proteomics :: differential-proteomics`
- `proteomics :: ms-proteomics-preprocessing`
- `proteomics :: sequence-to-function-modeling`
- `proteomics :: protein-embeddings`
- `transcriptomics :: spatial-transcriptomics`
- `transcriptomics :: pseudobulk-analysis`
- `transcriptomics :: rna-velocity`
- `genomics :: gwas`
- `genomics :: fine-mapping`
- `clinical_biomedical_data_science :: privacy-preserving-analysis`
- `systems_biology_and_network_science :: regulatory-network-inference`
- `proteomics :: protein-identification-quantification`

Validation for this phase:

- `python3 scripts/validate_repository.py`
- `python3 scripts/close_frontier_leaves.py`
- `python3 scripts/audit_skill_suite.py --json-out reports/refresh-20260321-cycle/audit_skill_suite.json --markdown-out reports/refresh-20260321-cycle/audit_skill_suite.md`
- `python3 scripts/run_skill_smoke_matrix.py --dry-run --json-out reports/refresh-20260321-cycle/smoke_matrix_dry_run.json --markdown-out reports/refresh-20260321-cycle/smoke_matrix_dry_run.md`
- `python3 scripts/benchmark_skill_advantage.py --case fastqc-multiqc-read-qc --json-out reports/refresh-20260321-cycle/benchmark_fastqc.json --markdown-out reports/refresh-20260321-cycle/benchmark_fastqc.md`
- `python3 scripts/build_site.py`

## Phase 54

Phase 54 was a conservative `tree_check` audit after the repository reached zero TODO and zero frontier leaves.

This pass did not add new resources or skills. Instead, it verified that the current post-closure tree remains internally consistent and that the next backlog is the unverified implementation queue.

Key findings:

- The tree stayed at `254` covered leaves, `0` frontier-only leaves, and `0` pure TODO leaves.
- The registry stayed at `366` canonical resources and `266` skills.
- Verified counts stayed at `116` `sandbox_verified` and `4` `slurm_verified`, while `146` skills still remain at `implemented`.
- `145` covered leaves still have zero verified skills, so tree-check success is now better measured by shrinking the unverified covered-leaf backlog than by frontier or TODO counts.
- `scripts/audit_skill_suite.py` remained clean on hard failures and `scripts/run_skill_smoke_matrix.py --dry-run` still resolved `266/266` smoke targets.
- A focused benchmark rerun for `fastqc-multiqc-read-qc` confirmed that the maintained wrapper still improves deliverable structure over the ad hoc baseline.

Refined next-pass queue for the requested priorities:

- `proteomics :: differential-proteomics`
- `proteomics :: ms-proteomics-preprocessing`
- `proteomics :: sequence-to-function-modeling`
- `proteomics :: protein-embeddings`
- `transcriptomics :: spatial-transcriptomics`
- `transcriptomics :: pseudobulk-analysis`
- `transcriptomics :: rna-velocity`
- `genomics :: gwas`
- `genomics :: fine-mapping`
- `clinical_biomedical_data_science :: privacy-preserving-analysis`
- `systems_biology_and_network_science :: regulatory-network-inference`
- `proteomics :: protein-identification-quantification`

The raw framework focus queue still behaves as designed, but broad terms such as `spatial` and `legal` are too coarse for this stage because they pull geospatial or other lexical matches ahead of spatial omics and privacy/regulatory leaves. The recommended follow-up is to use narrower focus terms during the next implementation cycle rather than changing framework behavior during this audit.

Validation for this phase:

- `python3 scripts/validate_repository.py`
- `python3 scripts/build_site.py`
- `python3 scripts/audit_skill_suite.py --json-out reports/tree-check-20260320-cycle/audit_skill_suite.json --markdown-out reports/tree-check-20260320-cycle/audit_skill_suite.md`
- `python3 scripts/run_skill_smoke_matrix.py --dry-run --json-out reports/tree-check-20260320-cycle/smoke_matrix_dry_run.json --markdown-out reports/tree-check-20260320-cycle/smoke_matrix_dry_run.md`
- `python3 scripts/sciskill_framework.py --json status --focus-limit 12 --focus-term protein --focus-term proteomics --focus-term spatial --focus-term sequence --focus-term law --focus-term legal --focus-term gwas --focus-term single-cell --focus-term "single cell"`
- `python3 scripts/benchmark_skill_advantage.py --case fastqc-multiqc-read-qc --json-out reports/tree-check-20260320-cycle/benchmark_fastqc.json --markdown-out reports/tree-check-20260320-cycle/benchmark_fastqc.md`

## Phase 43

Phase 43 completed the frontier-closure backlog by combining three paths:

- syncing previously implemented but unregistered verified skills,
- fixing one blocked runtime (`bcftools`) with a dedicated self-consistent prefix at `slurm/envs/bcftools`,
- generating deterministic implemented starters for every remaining frontier leaf so the full taxonomy is now skill-backed.

New verified-or-synced leaves included:

- `scientific_knowledge_access_and_method_extraction :: paper-triage-and-ranking`
- `scientific_knowledge_access_and_method_extraction :: review-paper-mining`
- `scientific_knowledge_access_and_method_extraction :: benchmark-paper-mining`
- `scientific_knowledge_access_and_method_extraction :: figure-table-extraction`
- `scientific_knowledge_access_and_method_extraction :: dataset-code-package-link-extraction`
- `data_acquisition_and_dataset_handling :: metadata-harmonization`
- `data_acquisition_and_dataset_handling :: synthetic-toy-dataset-generation-for-tests`
- `epigenomics_and_chromatin :: chromatin-interaction-hi-c`
- `genomics :: variant-filtering`

Generated implemented starters then closed every remaining frontier-only leaf, bringing the tree to:

- `254` covered leaves
- `0` frontier-only leaves
- `0` pure TODO leaves

Repository state after the full closure wave:

- `366` canonical resources
- `266` skills
- `116` `sandbox_verified`
- `4` `slurm_verified`
- `146` `implemented`

Validation for this phase:

- `make validate`
- `make build-site`
- `python3 -m unittest tests.smoke.test_phase43_frontier_completion_skills -v`
- `python3 -m unittest tests.integration.test_build_site -v`
- `python3 -m unittest tests.integration.test_skill_review_scripts -v`
- `make smoke-bcftools-filtering`
- `make smoke-semantic-scholar-triage`
- `make smoke-frontier-generated`
- `make smoke-deepchem-molgraph`
- `make test`

## Phase 42

Phase 42 was a resource-only deepening pass across already frontier-backed leaves in scientific knowledge extraction, genomics, epigenomics, proteomics, physics, agriculture, and HPC.

Added canonical resources:

- `openreview-py-github`
- `semantic-scholar-recommendations-api`
- `semantic-scholar-graph-api`
- `paperswithcode-data-github`
- `docling-github`
- `grobid-github`
- `kallisto-docs`
- `deepvariant-github`
- `bcftools-filtering-howto`
- `nf-core-sarek-usage`
- `cooler-docs`
- `meme-suite-overview-docs`
- `muon-docs`
- `openms-docs`
- `deepfri-github`
- `coffea-hep-docs`
- `jax-md-docs`
- `agml-docs`
- `alphasimr-introduction`
- `apsim-docs`
- `openmpi-docs`
- `dask-jobqueue-docs`
- `apptainer-bind-paths-docs`

Strengthened frontier leaves:

- `scientific_knowledge_access_and_method_extraction :: paper-triage-and-ranking`
- `scientific_knowledge_access_and_method_extraction :: review-paper-mining`
- `scientific_knowledge_access_and_method_extraction :: benchmark-paper-mining`
- `scientific_knowledge_access_and_method_extraction :: figure-table-extraction`
- `scientific_knowledge_access_and_method_extraction :: dataset-code-package-link-extraction`
- `genomics :: quantification`
- `genomics :: variant-calling`
- `genomics :: variant-filtering`
- `genomics :: somatic-pipelines`
- `epigenomics_and_chromatin :: chromatin-interaction-hi-c`
- `epigenomics_and_chromatin :: motif-analysis`
- `epigenomics_and_chromatin :: multiome-integration`
- `proteomics_and_protein_biology :: ms-proteomics-preprocessing`
- `proteomics_and_protein_biology :: sequence-to-function-modeling`
- `physics_and_astronomy :: high-energy-detector-data-analysis`
- `physics_and_astronomy :: differentiable-simulation`
- `agriculture_food_and_plant_science :: precision-agriculture-sensing`
- `agriculture_food_and_plant_science :: crop-genomics-and-breeding`
- `agriculture_food_and_plant_science :: yield-forecasting`
- `hpc_slurm_and_scaling :: multi-node-jobs`
- `hpc_slurm_and_scaling :: scratch-space-usage`

The tree stayed at `100 covered / 154 frontier / 0 TODO` leaves. That unchanged leaf-state is expected because the cycle improved resource depth for frontier leaves without adding new skills.

## Phase 41

Phase 41 converted nine additional frontier-only leaves into runnable verified starters using lightweight local package workflows.

Converted leaves:

- `data_acquisition_and_dataset_handling :: compression-decompression`
- `data_acquisition_and_dataset_handling :: format-conversion`
- `statistical_and_machine_learning_foundations_for_science :: statistical-testing`
- `statistical_and_machine_learning_foundations_for_science :: dimensionality-reduction`
- `statistical_and_machine_learning_foundations_for_science :: experimental-design`
- `computational_chemistry_and_molecular_simulation :: geometry-optimization`
- `meta_maintenance :: skill-deduplication`
- `meta_maintenance :: resource-deduplication`
- `visualization_and_reporting :: summary-pages-and-catalogs`

New verified skills:

- `numcodecs-compression-decompression-starter`
- `pyarrow-format-conversion-starter`
- `scipy-statistical-testing-starter`
- `umap-dimensionality-reduction-starter`
- `pydoe3-experimental-design-starter`
- `ase-geometry-optimization-starter`
- `rapidfuzz-skill-deduplication-starter`
- `datasketch-resource-deduplication-starter`
- `mkdocs-summary-catalog-starter`

This cycle also added three targeted canonical resources:

- `pyarrow-cookbook-docs`
- `rapidfuzz-docs`
- `datasketch-docs`

The tree moved from `91 covered / 163 frontier / 0 TODO` leaves to `100 covered / 154 frontier / 0 TODO` leaves.

## Phase 40

Phase 40 converted five frontier-only leaves into runnable verified starters and also fixed an alias-resolution bug that was undercounting existing systems-biology coverage in the generated tree.

Converted leaves:

- `scientific_knowledge_access_and_method_extraction :: citation-chaining`
- `systems_biology_and_network_science :: pathway-analysis`
- `systems_biology_and_network_science :: network-propagation`
- `computational_chemistry_and_molecular_simulation :: small-molecule-conformer-generation`
- `computational_chemistry_and_molecular_simulation :: force-field-assignment`

Alias-repaired leaves:

- `systems_biology_and_network_science :: reactome-identifier-enrichment`
- `systems_biology_and_network_science :: gene-set-tooling-from-bioconductor`

New verified skills:

- `openalex-citation-chain-starter`
- `reactome-pathway-analysis-starter`
- `networkx-network-propagation-starter`
- `rdkit-conformer-generation-starter`
- `openmm-forcefield-assignment-starter`

This cycle also added four targeted canonical resources:

- `openalex-work-object-docs`
- `openalex-filter-works-docs`
- `networkx-pagerank-docs`
- `openmm-forcefield-api-docs`

The tree moved from `84 covered / 170 frontier / 0 TODO` leaves to `91 covered / 163 frontier / 0 TODO` leaves.

## Phase 38

Phase 38 converted four frontier-only leaves into runnable verified starters with lightweight dedicated runtime reuse instead of new heavyweight toolchains.

Converted leaves:

- `scientific_computing_and_numerical_methods :: uncertainty-aware-simulation`
- `scientific_agents_and_automation :: evaluation-harnesses-for-scientific-agents`
- `robotics_lab_automation_and_scientific_instrumentation :: robotic-experiment-planning`
- `visualization_and_reporting :: dashboards`

New verified skills:

- `chaospy-uncertainty-propagation-starter`
- `inspect-evaluation-harness-starter`
- `autoprotocol-experiment-plan-starter`
- `dash-scientific-dashboard-starter`

This cycle also added four targeted canonical resources:

- `chaospy-github`
- `inspect-docs`
- `autoprotocol-python-github`
- `dash-basic-callbacks`

The tree moved from `80 covered / 174 frontier / 0 TODO` leaves to `84 covered / 170 frontier / 0 TODO` leaves.

## Phase 37

Phase 37 converted four frontier-only leaves into runnable verified starters while also repairing the responsive site layout.

Converted leaves:

- `data_acquisition_and_dataset_handling :: data-validation`
- `data_acquisition_and_dataset_handling :: data-provenance-tracking`
- `genomics :: alignment-and-mapping`
- `robotics_lab_automation_and_scientific_instrumentation :: instrument-control-and-scheduling`

New verified skills:

- `frictionless-tabular-validation-starter`
- `rocrate-metadata-bundle-starter`
- `minimap2-read-mapping-starter`
- `qcodes-parameter-sweep-starter`

This cycle also added three targeted canonical resources:

- `rocrate-py-github`
- `minimap2-manpage`
- `qcodes-measurement-tutorial`

The tree moved from `76 covered / 178 frontier / 0 TODO` leaves to `80 covered / 174 frontier / 0 TODO` leaves.

## Phase 35

Phase 35 converted five frontier-only leaves into runnable verified starters with minimal targeted runtime extension.

Converted leaves:

- `data_acquisition_and_dataset_handling :: chunking-sharding`
- `reproducible_workflows_and_workflow_engines :: reproducible-notebooks`
- `reproducible_workflows_and_workflow_engines :: ci-for-scientific-pipelines`
- `computational_chemistry_and_molecular_simulation :: molecular-dynamics-setup`
- `statistical_and_machine_learning_foundations_for_science :: bayesian-optimization`

New verified skills:

- `zarr-chunked-array-store-starter`
- `papermill-parameterized-notebook-starter`
- `github-actions-scientific-ci-starter`
- `openmm-langevin-dynamics-starter`
- `optuna-bayesian-optimization-starter`

This cycle also added three targeted canonical docs:

- `optuna-docs`
- `papermill-parameterize-docs`
- `github-actions-workflow-syntax-docs`

The tree moved from `71 covered / 183 frontier / 0 TODO` leaves to `76 covered / 178 frontier / 0 TODO` leaves.

## Phase 34

Phase 34 completed the resource-only TODO-elimination sweep.

- Added `101` new canonical resources without changing skill count.
- Eliminated all pure TODO leaves from `site/tree.json`.
- Shifted the tree from `71 covered / 82 frontier / 101 TODO` leaves to `71 covered / 183 frontier / 0 TODO` leaves.

This is the desired post-search state: the library now has resource coverage across the full explicit taxonomy, and the next work should focus on converting the highest-value frontier leaves into runnable skills.

## Phase 33

Phase 33 converted five resource-backed frontier leaves into runnable verified starters without adding new canonical resources.

Converted leaves:

- `genomics :: read-qc-and-trimming`
- `epigenomics_and_chromatin :: peak-calling`
- `proteomics_and_protein_biology :: sequence-feature-annotation`
- `neuroscience_and_neuroimaging :: connectomics-and-graph-analysis`
- `clinical_biomedical_data_science :: fairness-bias-analysis`

New verified skills:

- `fastqc-multiqc-read-qc-starter`
- `macs3-peak-calling-starter`
- `uniprot-sequence-feature-annotation-starter`
- `mne-connectivity-graph-starter`
- `fairlearn-bias-audit-starter`

The tree moved from `66 covered / 87 frontier / 101 TODO` leaves to `71 covered / 82 frontier / 101 TODO` leaves.
This cycle validates the new Phase 32 rule: once breadth seeding exposes many frontier leaves, the best next step is to convert the most executable leaves first.

## Phase 32

No top-level AI-for-Science domains are uncovered anymore. This cycle was intentionally resource-only and added `49` new canonical resources without changing skill count.

The tree moved from `66 covered / 0 frontier / 188 TODO` leaves to `66 covered / 87 frontier / 101 TODO` leaves. That change is desirable: many previously empty leaves now have canonical source material attached and are ready for implementation.

## Highest-Density Frontier Domains

| Area | Resources | Skills | Frontier leaves | Remaining TODO | Notes |
|------|-----------|--------|-----------------|----------------|-------|
| Scientific knowledge access | 13 | 5 | 7 | 3 | OpenReview, Papers with Code client, GROBID, Docling, and CODECHECK now seed triage, benchmark, protocol, reproducibility, and link-extraction leaves. |
| Data acquisition / dataset handling | 14 | 1 | 7 | 0 | The branch is now fully resource-backed through DataLad, DVC, Frictionless, Croissant, Arrow, Zarr, kerchunk, Great Expectations, RO-Crate, and numcodecs. |
| Genomics | 20 | 4 | 11 | 2 | The frontier now includes germline pipelines, GWAS, fine-mapping, PRS, comparative genomics, and phylogenomics via GATK, PLINK 2, susieR, PRSice, OrthoFinder, and IQ-TREE. |
| Epigenomics / chromatin | 14 | 1 | 9 | 0 | deepTools and SEACR now join MACS3, SnapATAC2, pycisTopic, Signac, ArchR, Bismark, TOBIAS, and HiCExplorer to seed nearly the whole branch. |
| Metabolomics / other omics | 10 | 1 | 5 | 1 | MetaboAnnotation, mixMC, and MetaboAnalyst now complement `xcms`, GNPS, QIIME 2, and mixOmics. |
| Proteomics / protein biology | 10 | 3 | 6 | 3 | UniProt, ESM, FragPipe, and PTM-Shepherd now push the branch past dataset lookup into annotation, embeddings, preprocessing, and PTM analysis. |
| Clinical / biomedical data science | 11 | 2 | 6 | 0 | OHDSI ATLAS and Fairlearn now deepen cohort extraction and fairness analysis alongside OMOP, MIMIC-IV, PyHealth, PyPOTS, AIF360, and Opacus. |
| Scientific agents / automation | 10 | 2 | 6 | 0 | OpenHands, Mem0, and paper-qa now extend code generation, long-horizon memory, and literature-to-hypothesis coverage beyond the existing orchestration and evaluation stack. |
| Physics / astronomy | 13 | 1 | 9 | 0 | Dedalus, NeuralOperator, SimPEG, PySR, Lightkurve, and Stingray now seed PDE, surrogate-model, inverse-problem, symbolic-discovery, and scientific time-series leaves. |
| Robotics / lab automation / instrumentation | 13 | 1 | 6 | 1 | Ophyd, Bluesky Queue Server, bluesky-adaptive, Databroker, pycro-manager, and BoTorch deepen control, LIMS-adjacent integration, microscopy, and closed-loop optimization. |
| Neuroscience / neuroimaging | 12 | 3 | 5 | 1 | SpikeInterface, MNE-Connectivity, DIPY, BrainGlobe Atlas API, and Nemos now seed spike sorting, connectomics, multimodal fusion, atlas integration, and neural decoding. |
| Scientific computing / numerical methods | 9 | 2 | 5 | 1 | dolfin-adjoint, Chaospy, Diffrax, and UQpy now seed inverse problems, uncertainty-aware simulation, and differentiable simulator programming. |

## Priority Queue

1. Convert frontier leaves with both multiple canonical resources and simple toy-data execution paths first: `bulk-rna-seq-qc-and-normalization`, `rna-velocity`, `pathway-analysis`, `virtual-screening`, `molecular-dynamics-setup`, `climate-reanalysis-access`, `bayesian-optimization`, and `ci-for-scientific-pipelines`.
2. Use the widest frontier domains as the next skill-conversion queue: `reproducible_workflows_and_workflow_engines`, `statistical_and_machine_learning_foundations_for_science`, `computational_chemistry_and_molecular_simulation`, `hpc_slurm_and_scaling`, `proteomics_and_protein_biology`, and `ecology_evolution_and_biodiversity`.
3. Run future resource-only passes only when new taxonomy branches are added or when a frontier leaf still lacks a sufficiently canonical anchor.

## 2026-03-21 Refresh Note

- Refresh re-ran the shared validator, frontier-closure, audit, dry-run smoke-matrix, and focused benchmark surfaces for the active cycle.
- Repository state remained stable at `366` deduplicated resources, `266` skills, `254` covered leaves, `0` frontier leaves, and `0` TODO leaves.
- `scripts/audit_skill_suite.py` still reports `266/266` smoke-covered skills, `253/266` skills with local tests, `116` `sandbox_verified` skills, `4` `slurm_verified` skills, and `0` hard failures.
- `scripts/run_skill_smoke_matrix.py --dry-run` still resolves smoke targets for all `266` registered skills with no missing targets.
- The focused refresh benchmark kept the maintained wrappers ahead on deliverable rate for both GWAS and interface-analysis starter cases while command success stayed tied at `1.0`.
- Refresh completed the missing `manifest.json` and `stages/05-refresh/result.json` for `reports/framework-runs/20260321-042804-protein-spatial-sequence-law-gwas-single-cell-cycle-006`, allowing site generation to index that cycle conservatively without changing verification labels.

## 2026-03-21 Tree Check Note

- Tree check re-ran the shared validator, site build, frontier-closure sync, skill-suite audit, dry-run smoke matrix, and focused GWAS/MS-proteomics benchmark surfaces.
- Repository state remained stable at `366` deduplicated resources, `266` skills, `254` covered leaves, `0` frontier leaves, and `0` TODO leaves.
- `scripts/close_frontier_leaves.py` remained a no-op, confirming the tree is fully closed under the current taxonomy.
- `scripts/audit_skill_suite.py` still reports `0` hard failures, `266/266` smoke-covered skills, and `253/266` skills with local tests.
- `scripts/run_skill_smoke_matrix.py --dry-run` still resolves smoke targets for all `266` registered skills with no missing mappings.
- The focused benchmark pass kept maintained-wrapper deliverable rate at `1.0` for both GWAS and MS proteomics preprocessing starter cases, versus `0.625` aggregate for the ad hoc baseline.
- The refined next-pass queue is now centered on unverified leaves in proteomics, spatial omics, GWAS, and privacy/regulatory follow-on work rather than already-verified single-cell preprocessing and integration.
