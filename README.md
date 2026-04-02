# SciSkillUniverse

SciSkillUniverse is a structured repository for reusable AI-for-science skills spanning computational biology, proteomics, cheminformatics, workflow automation, and adjacent computational-science tooling. The repository is organized around a registry, skill folders with provenance and tests, lightweight site artifacts, and infrastructure notes for workflow engines and Slurm.

## Repository Layout

- `registry/`: taxonomy, resource registry, skill index, coverage, and update logs
- `skills/`: reusable skill folders grouped by domain
- `tests/`: smoke, regression, integration, and Slurm-oriented checks
- `site/`: generated JSON for browsing plus a static `index.html`
- `slurm/`: environment, job, log, and report scaffolding
- `scripts/`: dependency-light validators and site generators
- `reports/`: per-run summaries

## Quickstart

```bash
make validate
make build-site
make test
```

Additional smoke checks:

```bash
make smoke-openalex
make smoke-europepmc
make smoke-encode
make smoke-biosamples
make smoke-metabolights
make smoke-quickgo
make smoke-ncbi-gene
make smoke-literature-brief
make smoke-reactome
make smoke-reactome-enrichment
make smoke-pride
make smoke-ebi-proteins
make smoke-chembl
make smoke-cwl
make smoke-wdl
make smoke-rdkit
make smoke-rdkit-standardize
make smoke-deepchem
make smoke-openmm
make smoke-psi4
make smoke-fgsea
make smoke-clusterprofiler
make smoke-snakemake
make smoke-nextflow
make smoke-nextflow-slurm
make smoke-nf-core
make smoke-scanpy
make smoke-scanpy-ranked-genes
make smoke-pydeseq2
make smoke-sourmash
make smoke-skimage
make smoke-skimage-regionprops
make smoke-xarray
make smoke-rasterio
make smoke-matminer
make smoke-matminer-regression
make smoke-geopandas
make smoke-pymc
make smoke-arviz
make smoke-dowhy
make smoke-skill-router
make smoke-matplotlib
make smoke-plotly
make smoke-quarto
make smoke-link-audit
make smoke-precommit
make smoke-slurm
```

## Automation Framework

The repository includes a `codex exec` orchestration layer for the recurring `tree_check -> resource_search -> skill_build -> skill_test -> refresh` loop. The CLI lives in [`scripts/sciskill_framework.py`](/home/shuaikes/projects/agent/SciSkillUniverse/scripts/sciskill_framework.py) and reuses the existing registry, site builder, validation, smoke, audit, and benchmark scripts instead of introducing a second control plane.

### Commands

Use `status` to inspect the current repository summary and the next leaves worth working on:

```bash
python3 scripts/sciskill_framework.py --json status --focus-limit 10
```

Use `cycle` to run one or more automation loops. By default it executes the staged flow defined in the framework prompts:

```bash
python3 scripts/sciskill_framework.py cycle --loops 1 --verification-mode standard
python3 scripts/sciskill_framework.py cycle --loops 2 --focus-limit 12 --label nightly-pass
python3 scripts/sciskill_framework.py cycle \
  --stages tree_check,resource_search,skill_build,skill_test,refresh \
  --extra-context "Prioritize uncovered leaves in robotics and physics."
```

Use `design-skill` when you already have a specific task or user prompt and want the framework to design a skill around it:

```bash
python3 scripts/sciskill_framework.py design-skill \
  --prompt "Design a skill for literature-backed pathway enrichment benchmarking." \
  --verification-mode validate
```

Use `evaluate-skills` to run the hierarchical evaluation framework against existing skills, including correctness repair, benchmark comparison against a no-skill baseline, and novelty checking:

```bash
python3 scripts/sciskill_framework.py evaluate-skills \
  --skill-slug openalex-literature-search \
  --layer1-fix-attempts 2 \
  --layer2-optimize-attempts 2 \
  --verification-mode validate

python3 scripts/sciskill_framework.py evaluate-skills \
  --all \
  --verification-mode none \
  --label full-library-review
```

Recent layer-2 benchmarks live under `reports/framework-runs/20260321-031855-protein-spatial-sequence-law-gwas-single-cell-eval-004/skills/gbif-species-occurrence-search-starter/layer2/benchmark-00/` for GBIF species occurrence search and `reports/framework-runs/20260321-072208-protein-spatial-sequence-law-gwas-single-cell-eval-010/skills/semantic-scholar-review-paper-mining-starter/layer2/benchmark-00/` for Semantic Scholar review-paper mining. Both compare the maintained wrapper against a minimal ad hoc baseline on concrete task-specific cases.
The latest Time-series clinical modeling Starter layer-2 benchmark lives under `reports/framework-runs/20260321-080211-protein-spatial-sequence-law-gwas-single-cell-eval-011/skills/time-series-clinical-modeling-starter/layer2/benchmark-00/` and compares the maintained starter wrapper against a note-only baseline on canonical-summary, augmented-objective, and nested-output cases.
The latest I/O-aware workflows layer-2 benchmark lives under `reports/framework-runs/20260321-031855-protein-spatial-sequence-law-gwas-single-cell-eval-004/skills/i-o-aware-workflows-starter/layer2/benchmark-00/` and compares the maintained starter wrapper against a no-skill summary writer on three deterministic local cases.
The latest NCBI Gene Search layer-2 benchmark lives under `reports/framework-runs/20260321-044031-protein-spatial-sequence-law-gwas-single-cell-eval-006/skills/ncbi-gene-search/layer2/benchmark-00/` and compares the maintained wrapper against a thin ad hoc projection built from the same local fixture payloads.
The latest OpenAlex Citation Chain layer-2 benchmark lives under `reports/framework-runs/20260321-051402-protein-spatial-sequence-law-gwas-single-cell-eval-007/skills/openalex-citation-chain-starter/layer2/benchmark-00/` and compares the maintained wrapper against an ad hoc upstream-only citation summary on DOI-based citation-chain cases.
The latest PRIDE Project Search layer-2 benchmark lives under `reports/framework-runs/20260321-051402-protein-spatial-sequence-law-gwas-single-cell-eval-007/skills/pride-project-search/layer2/benchmark-00/` and compares the maintained wrapper against a compact ad hoc JSON writer on canonical, noisy-normalization, and sparse-fallback PRIDE snapshot cases.
The latest Population dynamics and ecological forecasting starter layer-2 benchmark lives under `reports/framework-runs/20260321-051402-protein-spatial-sequence-law-gwas-single-cell-eval-007/skills/population-dynamics-and-ecological-forecasting-starter/layer2/benchmark-00/` and compares the maintained starter wrapper against an ad hoc notes writer on canonical and augmented local starter-plan cases.
The latest OpenMM System Minimization layer-2 benchmark lives under `reports/framework-runs/20260321-051402-protein-spatial-sequence-law-gwas-single-cell-eval-007/skills/openmm-system-minimization/layer2/benchmark-00/` and compares the maintained wrapper against a direct ad hoc OpenMM minimizer on canonical completeness and nested-output-path robustness cases.
The latest Perturb-seq Starter layer-2 benchmark lives under `reports/framework-runs/20260321-051402-protein-spatial-sequence-law-gwas-single-cell-eval-007/skills/perturb-seq-starter/layer2/benchmark-00/` and compares the maintained starter wrapper against a minimal ad hoc notes writer on canonical and augmented local starter-plan cases.
The latest Phylogenomics Starter layer-2 benchmark lives under `reports/framework-runs/20260321-051402-protein-spatial-sequence-law-gwas-single-cell-eval-007/skills/phylogenomics-starter/layer2/benchmark-00/` and compares the maintained starter wrapper against a minimal ad hoc notes writer on canonical, augmented, and nested-output starter-plan cases.
The latest Operator Learning and Surrogate Physics Models starter benchmark lives under `reports/framework-runs/20260321-051402-protein-spatial-sequence-law-gwas-single-cell-eval-007/skills/operator-learning-and-surrogate-physics-models-starter/layer2/benchmark-00/` and compares the maintained starter wrapper against a registry-derived ad hoc note writer on nested-contract, promotion-guidance, and handoff-completeness cases.
The latest PlantCV plant-phenotyping layer-2 benchmark lives under `reports/framework-runs/20260321-051402-protein-spatial-sequence-law-gwas-single-cell-eval-007/skills/plantcv-plant-phenotyping-starter/layer2/benchmark-00/` and compares the maintained wrapper against a minimal NumPy/imageio baseline on canonical, nested-output, and threshold-propagation toy-plant cases.
The latest Precision agriculture sensing Starter layer-2 benchmark lives under `reports/framework-runs/20260321-051402-protein-spatial-sequence-law-gwas-single-cell-eval-007/skills/precision-agriculture-sensing-starter/layer2/benchmark-00/` and compares the maintained wrapper against a minimal ad hoc notes writer on canonical, nested-output, and augmented starter-plan cases.
The latest Reactome Pathway Analysis Starter layer-2 benchmark lives under `reports/framework-runs/20260321-063711-protein-spatial-sequence-law-gwas-single-cell-eval-009/skills/reactome-pathway-analysis-starter/layer2/benchmark-00/` and compares the maintained wrapper against a no-skill direct Reactome fetch baseline on canonical, noisy-inline, and noisy-file-input cases.
The latest RO-Crate Metadata Bundle Starter layer-2 benchmark lives under `reports/framework-runs/20260321-063711-protein-spatial-sequence-law-gwas-single-cell-eval-009/skills/rocrate-metadata-bundle-starter/layer2/benchmark-00/` and compares the maintained wrapper against a minimal ad hoc RO-Crate build baseline on canonical, custom-metadata, and stale nested-output cases.
The latest Scanpy QC Starter layer-2 benchmark lives under `reports/framework-runs/20260321-063711-protein-spatial-sequence-law-gwas-single-cell-eval-009/skills/scanpy-qc-starter/layer2/benchmark-00/` and compares the maintained wrapper against a raw tabular summary baseline on canonical and augmented mitochondrial-rich toy count matrices.
The latest Sparse / iterative linear algebra Starter layer-2 benchmark lives under `reports/framework-runs/20260321-072208-protein-spatial-sequence-law-gwas-single-cell-eval-010/skills/sparse-iterative-linear-algebra-starter/layer2/benchmark-00/` and compares the maintained wrapper against a minimal starter-note baseline on canonical, augmented-objective, and nested-output cases.
The latest Spatial interpolation and uncertainty Starter layer-2 benchmark lives under `reports/framework-runs/20260321-072208-protein-spatial-sequence-law-gwas-single-cell-eval-010/skills/spatial-interpolation-and-uncertainty-starter/layer2/benchmark-00/` and compares the maintained wrapper against a notes-only baseline on canonical, checklist-propagation, and nested-output cases.

Use `campaign` for a long-running checkpointable job that alternates targeted discovery/build cycles with batched whole-library evaluation:

```bash
python3 scripts/sciskill_framework.py \
  --workspace-root scratch/framework/workspaces-$SLURM_JOB_ID \
  campaign \
  --label protein-spatial-sequence-law-gwas-single-cell \
  --focus-term protein \
  --focus-term proteomics \
  --focus-term spatial \
  --focus-term sequence \
  --focus-term law \
  --focus-term legal \
  --focus-term gwas \
  --focus-term single-cell \
  --focus-term "single cell" \
  --focus-limit 12 \
  --stage-workers 6 \
  --background-validation-limit 24 \
  --background-validation-workers 6 \
  --evaluation-batch-size 24 \
  --evaluation-workers 6 \
  --verification-mode standard \
  --max-runtime-minutes 450 \
  --stop-buffer-minutes 20
```

Use `campaign-status` to inspect the latest checkpoint summary without reopening the full manifests:

```bash
python3 scripts/sciskill_framework.py campaign-status \
  --label protein-spatial-sequence-law-gwas-single-cell
```

Interactive runs now emit live progress to the terminal on `stderr`, including stage banners, current skill / command, command tails, and heartbeat lines while a long `codex exec` stage is still running. This is intentionally separate from the final payload on `stdout`.
`[STDERR]` now means a subprocess wrote to the `stderr` stream; it does not automatically mean the stage failed. Real framework failures are logged as `[ERR ]`, and soft fallbacks or retries stay at `[WARN]`.

Global options such as `--codex-bin`, `--model`, `--reasoning-effort`, `--profile`, `--full-auto`, `--stage-timeout`, `--verification-timeout`, and `--workspace-root` apply to all subcommands. This makes it easy to point the framework at a different Codex binary, a different fallback model/profile, or a CPU-node-local worker workspace without editing repo code.

### Default Stage Model Routing

By default, the framework now routes different phases to different `codex exec` model / effort pairs to control token spend:

- `tree_check` -> `gpt-5.4` + `medium`
- `resource_search` -> `gpt-5.4` + `high`
- `skill_build` -> `gpt-5.4` + `medium`
- `skill_test` -> `gpt-5.4-mini` + `medium`
- `design_skill` -> `gpt-5.4` + `medium`
- `layer1_fix`, `layer2_benchmark`, `layer2_optimize`, `novelty_check` -> `gpt-5.4-mini` + `medium`
- `refresh` -> `gpt-5.4` + `medium`

The global `--model` / `--reasoning-effort` flags still work as fallbacks for stages that do not have an explicit route.

### Parallel Execution

The framework now supports two independent forms of parallelism:

- leaf-parallel discovery/build/test inside `cycle` with `--stage-workers N`
- per-layer parallel skill evaluation inside `evaluate-skills` with `--workers N`

`resource_search`, `skill_build`, and `skill_test` can fan out into multiple isolated worker workspaces. Those workers run `codex exec` concurrently, keep their writes local during the stage, and only sync back leaf-owned paths such as `skills/**`, leaf-specific `tests/test_*.py`, and `slurm/jobs/**` before `refresh`. Shared files like `registry/**`, `site/**`, `README.md`, `experiments.md`, and `Makefile` remain sequential and are updated in the final refresh path.

`cycle` can also run a background validation lane at the same time as discovery/build workers:

```bash
python3 scripts/sciskill_framework.py --workspace-root scratch/framework/workspaces-$SLURM_JOB_ID cycle \
  --loops 1 \
  --focus-limit 12 \
  --stage-workers 4 \
  --background-validation-limit 16 \
  --background-validation-workers 4 \
  --verification-mode standard \
  --label parallel-cycle
```

That command runs concurrent Codex workers for leaf stages while separately validating existing skills in the main repository. It is the recommended way to let discovery/development and validation progress simultaneously without racing on shared files.
Use that single `cycle` run for concurrent discovery plus validation; launching multiple top-level framework commands at the same time is still discouraged because the live dashboard intentionally tracks one active run at a time.

### Typical Workflow

1. Inspect the current focus queue with `status`.
2. Run a bounded cycle with `cycle --loops 1`.
3. Review the generated manifest and stage artifacts under `reports/framework-runs/`.
4. Rebuild the site with `make build-site` if you want the framework dashboard to reflect the new run immediately.
5. Open [`site/index.html`](/home/shuaikes/projects/agent/SciSkillUniverse/site/index.html) to inspect the updated run deck, skill catalog, and taxonomy state.

The default stage semantics are:

- `tree_check`: inspect taxonomy coverage and identify the next high-value leaves.
- `resource_search`: search for papers, docs, repos, notebooks, workflows, and API references before implementation.
- `skill_build`: construct or extend skills from the gathered evidence.
- `skill_test`: add or repair verification, smoke checks, and Slurm-backed tests when appropriate.
- `refresh`: update registry artifacts, site data, reports, and planning files.

### Hierarchical Evaluation

`evaluate-skills` is the standalone entry point for testing newly built skills and re-reviewing the historical skill library. It implements three automated layers and leaves the human layer manual:

1. `layer1 correctness`: run the resolved smoke or registry test command for each selected skill. If the check fails, the framework calls `codex exec` in a `layer1_fix` loop, captures debug findings, applies fixes, and reruns the check.
2. `layer2 benchmark`: for skills that pass layer 1, the framework calls `codex exec` to design task-specific benchmark cases and compare a with-skill path against a no-skill baseline. If the skill does not beat the baseline clearly enough, the framework enters a `layer2_optimize` loop, reruns correctness, and benchmarks again.
3. `novelty check`: the framework asks `codex exec` to compare the skill against local similarity candidates and external tool ecosystems, explicitly including ToolUniverse overlap review.

Layer-2 benchmark artifacts are written under `reports/framework-runs/.../skills/<skill-slug>/layer2/benchmark-00/` and include the JSON manifest plus a short Markdown summary for the selected cases.
For starter-style skills, the benchmark should emphasize the actionability fields added by the maintained wrapper, especially the starter steps and promotion checklist.

Layer-2 benchmark runs leave their case notes and JSON outputs under `reports/framework-runs/<run-id>/skills/<skill-slug>/layer2/`, so the comparison can be audited after the fact.

Human expert evaluation is intentionally not automated. The manifest records this as a manual layer so the hierarchy stays explicit.

Parallel evaluation is enabled with `--workers`:

```bash
python3 scripts/sciskill_framework.py evaluate-skills \
  --skill-slug openalex-literature-search \
  --skill-slug europepmc-method-triage \
  --workers 4 \
  --verification-mode validate \
  --label parallel-eval
```

`evaluate-skills` parallelizes layer 1 correctness, layer 2 benchmark, and novelty review across the selected skills. It does not run those three layers simultaneously with each other; each layer stays ordered, but each layer can fan out to multiple skill workers.

`campaign` adds checkpointable batch orchestration on top of those primitives. Each campaign iteration runs:

1. one targeted `cycle` with parallel discovery/build/test workers;
2. one batched `evaluate-skills` pass over the next pending skills that have not already been evaluated in that campaign state.

Completed evaluation results are written into the campaign checkpoint and are skipped on resume. In addition, interrupted evaluation batches are now resumed from existing per-skill `layer1_record.json`, `layer2_record.json`, and `novelty_record.json` artifacts, so a killed CPU-node job does not have to restart the whole batch from the beginning.

### Verification Modes

Post-stage verification is explicit and conservative:

- `none`: skip post-run verification entirely.
- `validate`: run [`scripts/validate_repository.py`](/home/shuaikes/projects/agent/SciSkillUniverse/scripts/validate_repository.py).
- `standard`: run validation plus [`scripts/build_site.py`](/home/shuaikes/projects/agent/SciSkillUniverse/scripts/build_site.py).
- `full`: run validation, site generation, and the full `unittest` suite.
- `audit`: extend `full` with [`scripts/audit_skill_suite.py`](/home/shuaikes/projects/agent/SciSkillUniverse/scripts/audit_skill_suite.py) and a dry-run smoke-matrix pass.

Use `validate` or `standard` while iterating quickly, and reserve `full` or `audit` for larger refactors or batch skill additions.

### Artifacts and Website Integration

The framework writes durable run artifacts to:

- `reports/framework-runs/<timestamp>-<label>/manifest.json`: top-level record of the run, requested stages, verification mode, and before/after summaries.
- `reports/framework-runs/<timestamp>-<label>/stages/*/`: per-stage Codex prompts, outputs, and normalized result records.
- `scratch/framework/latest_status.json`: latest repository snapshot from `status`.
- `scratch/framework/latest_active_run.json`: live run-state snapshot used by the website while a cycle, design pass, or evaluation is still executing.
- `scratch/framework/latest_run.json` and `scratch/framework/latest_design_skill.json`: latest run-level state for the most recent automation actions.
- `scratch/framework/latest_evaluation.json`: latest hierarchical evaluation manifest for the standalone evaluation entrypoint.
- `scratch/framework/campaigns/<label>/campaign_state.json`: checkpointable campaign state with completed evaluation slugs, cycle runs, and aggregated counters.
- `scratch/framework/campaigns/<label>/status.json`: compact machine-readable status for monitoring active or resumed campaigns.
- `scratch/framework/campaigns/<label>/summary.md`: human-readable aggregate table with designed skills, new skills, layer-1 results, and layer-2 results.
- `scratch/framework/campaigns/<label>/campaign.log`: durable outer-loop campaign log.

Hierarchical evaluation runs also write per-skill artifacts under `reports/framework-runs/<timestamp>-<label>/skills/<skill-slug>/`, including layer-1 checks, fix attempts, benchmark attempts, optimization feedback, and novelty summaries.

The site generator now also emits [`site/framework_runs.json`](/home/shuaikes/projects/agent/SciSkillUniverse/site/framework_runs.json). The framework dashboard inside [`site/index.html`](/home/shuaikes/projects/agent/SciSkillUniverse/site/index.html) reads that file to show the latest run, stage ledger, and focus queue beside the skill browser. A typical refresh loop is:

```bash
python3 scripts/sciskill_framework.py cycle --loops 1 --verification-mode standard
make build-site
python3 -m http.server 8000
```

Then open `http://localhost:8000/site/` and jump to the framework run deck.

While a run is still active, the website also polls `scratch/framework/latest_active_run.json` and renders a live shell with:

- a running / completed / failed status badge
- a spinner beside the currently executing stage
- stage chips for discovery, build, test, benchmark, and verification progress
- color-coded terminal lines for commands, stdout, neutral `stderr`, warnings, real errors, and success messages

Campaign monitoring is file-based rather than website-native today. `campaign-status` now reconciles existing `cycle-*` / `eval-*` artifacts before printing, so it reflects interrupted-but-recoverable runs more accurately. The recommended live views are:

```bash
python3 scripts/sciskill_framework.py campaign-status --label protein-spatial-sequence-law-gwas-single-cell
tail -f scratch/framework/campaigns/protein-spatial-sequence-law-gwas-single-cell/campaign.log
```

### CPU-Node Operation

For single-worker inspection, the login/head node is usually sufficient. For parallel discovery or parallel evaluation, prefer a CPU allocation so the control plane, `codex exec` processes, smoke commands, and worker workspace copies do not overload the login node.

Interactive example:

```bash
salloc -p cpu -N 1 -n 1 --cpus-per-task=8 --mem=32G -t 08:00:00
codex login status
python3 scripts/sciskill_framework.py --workspace-root scratch/framework/workspaces-$SLURM_JOB_ID cycle \
  --loops 1 \
  --focus-limit 12 \
  --stage-workers 4 \
  --background-validation-limit 16 \
  --background-validation-workers 4 \
  --verification-mode standard \
  --label cpu-parallel-cycle
```

If `codex login status` is not authenticated on the allocated node, run `codex login` there before starting the framework. Slurm itself does not need a separate Codex authorization step; the important requirement is that the node executing `scripts/sciskill_framework.py` has a valid Codex login and the same repo checkout.

For the long targeted campaign in this repository, the provided batch entrypoint is:

```bash
sbatch slurm/jobs/framework_domain_campaign_cpu.sbatch
```

That batch script already does `nvm use 24`, keeps CPU usage modest (`8` CPUs), writes Slurm logs under `slurm/logs/`, and uses the checkpointable `campaign` entrypoint so walltime-limited runs can be resumed safely. If a batch job is stopped mid-evaluation, resubmitting the same label resumes from the existing evaluation artifacts instead of re-running finished layer-1 / layer-2 work.

### Full-Repository Launch Examples

For a whole-repository discovery pass with standard verification:

```bash
python3 scripts/sciskill_framework.py --workspace-root scratch/framework/workspaces-$SLURM_JOB_ID cycle \
  --loops 1 \
  --focus-limit 254 \
  --stage-workers 8 \
  --background-validation-limit 32 \
  --background-validation-workers 8 \
  --verification-mode standard \
  --label full-repo-discovery
```

For a whole-library hierarchical evaluation pass over the existing skills:

```bash
python3 scripts/sciskill_framework.py evaluate-skills \
  --all \
  --workers 8 \
  --verification-mode validate \
  --layer1-fix-attempts 1 \
  --layer2-optimize-attempts 1 \
  --label full-library-eval
```

`evaluate-skills --all` is a long-running batch job. With the current library size, it can take hours because correctness, benchmark, and novelty review may each touch every selected skill. Use a smaller scoped pass first when debugging:

```bash
python3 scripts/sciskill_framework.py evaluate-skills \
  --skill-slug openalex-literature-search \
  --verification-mode validate \
  --layer1-fix-attempts 0 \
  --layer2-optimize-attempts 0 \
  --label smoke-eval-openalex
```

If you prefer the shorter Make wrappers, use:

```bash
make framework-status
make framework-cycle
make framework-cycle-parallel
make framework-evaluate-skills
make framework-evaluate-skills-parallel
make framework-campaign
make framework-campaign-status
make framework-submit-campaign
```

## Current Seeded Skills

- `openalex-literature-search`: live literature search via the OpenAlex API
- `europepmc-method-triage`: life-science literature triage via Europe PMC
- `crossref-metadata-search`: DOI-oriented metadata search via Crossref
- `ncbi-pubmed-search`: PubMed ID search via NCBI E-utilities
- `multi-source-literature-brief`: normalized cross-source brief across OpenAlex, Europe PMC, Crossref, and PubMed
- `ensembl-gene-lookup`: gene symbol lookup via the Ensembl REST API
- `ncbi-gene-search`: Entrez Gene symbol search and compact NCBI gene summaries
- `rcsb-pdb-search`: protein structure ID search via the RCSB PDB Search API
- `rcsb-pdb-entry-summary`: metadata lookup for a known PDB accession via the RCSB data API
- `clinicaltrials-study-search`: live trial registry search via ClinicalTrials.gov API v2
- `reactome-event-summary`: stable-ID summary lookup through the Reactome Content Service
- `reactome-identifiers-enrichment`: compact identifier-list pathway enrichment through the Reactome Analysis Service
- `pride-project-search`: PRIDE Archive project discovery for proteomics datasets
- `ebi-proteins-entry-summary`: accession-based protein summary lookup through the EBI Proteins API
- `chembl-molecule-search`: compact compound lookup through ChEMBL web services
- `encode-experiment-search`: public ENCODE experiment search for assay and biosample metadata
- `biosamples-sample-search`: public EBI BioSamples metadata search by free text
- `metabolights-study-search`: MetaboLights study discovery with accession-level summaries
- `quickgo-term-search`: GO term search through the official QuickGO ontology API
- `cwl-commandlinetool-starter`: runnable CWL CommandLineTool starter via the repo-managed `cwltool` prefix
- `wdl-task-starter`: runnable WDL greeting workflow via the repo-managed `miniwdl` + `udocker` prefix
- `rdkit-molecular-descriptors`: local RDKit descriptor summary generation from one SMILES string
- `rdkit-molecule-standardization`: local RDKit MolStandardize cleanup, fragment-parent, and canonical-tautomer summary generation
- `deepchem-circular-featurization`: local DeepChem circular fingerprint generation through a dedicated Python 3.9 prefix
- `openmm-system-minimization`: deterministic OpenMM toy-system minimization with energy-drop summary
- `psi4-single-point-energy`: deterministic Psi4 water single-point energy through the dedicated quantum-chemistry prefix
- `fgsea-preranked-enrichment`: local pre-ranked enrichment analysis via the repo-managed Bioconductor prefix
- `clusterprofiler-custom-enrichment`: local over-representation enrichment with custom TERM2GENE/TERM2NAME tables
- `snakemake-toy-workflow-starter`: runnable toy Snakemake workflow with deterministic outputs
- `nextflow-hello-workflow`: self-contained Nextflow verification workflow for both local and Slurm executor paths
- `nf-core-pipeline-list`: cleaned machine-readable nf-core catalog listing
- `slurm-job-debug-template`: reusable Slurm smoke-job patterns with real cluster accounting capture
- `scanpy-qc-starter`: toy Scanpy QC execution with AnnData export
- `scanpy-ranked-genes-starter`: toy Scanpy marker ranking across grouped cells
- `scanpy-dpt-trajectory-starter`: deterministic Scanpy diffusion-pseudotime trajectory starter
- `pydeseq2-differential-expression-starter`: deterministic bulk RNA-seq differential-expression starter with PyDESeq2
- `sourmash-signature-compare-starter`: deterministic MinHash sketch comparison for two toy DNA sequences
- `skimage-otsu-segmentation-starter`: deterministic toy segmentation with component summaries
- `skimage-regionprops-feature-extraction`: deterministic morphology and intensity feature extraction via `regionprops_table`
- `xarray-climate-cube-starter`: deterministic climate-cube summarization via Xarray
- `rasterio-windowed-raster-preprocessing-starter`: deterministic window extraction and average-resampling starter for toy raster preprocessing
- `matminer-composition-featurization`: stoichiometry-style composition features for simple formulas
- `matminer-property-regression-starter`: toy materials-property regression built from matminer composition features and Ridge regression
- `geopandas-spatial-join-starter`: point-in-polygon starter with built-in `PROJ` / `GDAL` path repair
- `pymc-bayesian-linear-regression-starter`: Bayesian linear regression with posterior means and 90% intervals
- `arviz-posterior-diagnostics-starter`: deterministic posterior-diagnostics summary with R-hat and effective sample sizes
- `dowhy-average-treatment-effect-starter`: deterministic backdoor ATE estimation plus placebo refutation with DoWhy
- `skill-registry-router-starter`: deterministic free-text routing from task queries to candidate skills in the local registry
- `matplotlib-publication-plot-starter`: deterministic two-panel Matplotlib figure rendering plus compact fit summary
- `plotly-interactive-report-starter`: self-contained interactive Plotly HTML report plus compact trend summary
- `quarto-notebook-report-starter`: deterministic notebook-to-HTML conversion via a repo-managed Quarto wrapper
- `registry-link-audit-starter`: dependency-light reachability audit for selected registry resources
- `precommit-regression-testing-starter`: deterministic tiny-repository `pre-commit` regression harness for local hook execution

## Design Notes

- Registry records are stored as JSONL for deterministic validation and easy diffs.
- `metadata.yaml` files use YAML-compatible JSON so the repository stays dependency-light.
- Heavy toolchains live in dedicated prefixes under `slurm/envs/` so Java/Nextflow/nf-core, Scanpy, workflow-language runtimes, chemistry toolkits, DeepChem, Psi4, Bioconductor packages, GeoPandas, and PyMC do not contaminate the base shell.
- The site browser now renders the full 21-domain taxonomy, and all 21 top-level domains now have at least one registered resource or skill.
- Verification labels are conservative but now broader: `50` skills are `sandbox_verified`, `2` are `slurm_verified`, and the current registry spans `52` skills and `70` canonical resources.

The latest `interactive-debug-workflows-starter` layer-2 benchmark lives under `reports/framework-runs/20260321-031855-protein-spatial-sequence-law-gwas-single-cell-eval-004/skills/interactive-debug-workflows-starter/layer2/benchmark-00/` and compares the maintained starter wrapper against a no-skill summary writer on three deterministic local cases.
The latest `plotly-interactive-report-starter` layer-2 benchmark lives under `reports/framework-runs/20260321-051402-protein-spatial-sequence-law-gwas-single-cell-eval-007/skills/plotly-interactive-report-starter/layer2/benchmark-00/fresh-run/` and compares the maintained wrapper against an ad hoc Plotly writer on canonical, flat-line, and nested-output cases.
The latest `pydoe3-experimental-design-starter` layer-2 benchmark lives under `reports/framework-runs/20260321-055525-protein-spatial-sequence-law-gwas-single-cell-eval-008/skills/pydoe3-experimental-design-starter/layer2/benchmark-00/` and compares the maintained wrapper against a compact ad hoc design writer on canonical two-factor and augmented three-factor cases.
The latest `pydeseq2-differential-expression-starter` layer-2 benchmark lives under `reports/framework-runs/20260321-055525-protein-spatial-sequence-law-gwas-single-cell-eval-008/skills/pydeseq2-differential-expression-starter/layer2/benchmark-00/` and compares the maintained wrapper against a raw fold-change baseline on canonical, reversed-contrast, and shuffled-metadata cases.
The latest `scanpy-ranked-genes-starter` layer-2 benchmark lives under `reports/framework-runs/20260321-063711-protein-spatial-sequence-law-gwas-single-cell-eval-009/skills/scanpy-ranked-genes-starter/layer2/benchmark-00/` and compares the maintained wrapper against an ad hoc Scanpy marker-table writer on canonical, deeper-summary, and missing-group-label cases.
