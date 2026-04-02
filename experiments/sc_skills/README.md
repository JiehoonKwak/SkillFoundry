# sc_skills Experiments

This directory is the experiment-only workspace for portable single-cell and spatial skills derived from the Spatial Agent task references in `ref/skill/`.

## Scope

- Source references live in [`ref/skill/`](/home/shuaikes/projects/agent/SciSkillUniverse/ref/skill).
- Generated experiment skills stay under [`experiments/sc_skills/`](/home/shuaikes/projects/agent/SciSkillUniverse/experiments/sc_skills).
- These experiment skills are not production skills and are not registered in `registry/`.
- The batch design flow can optionally route a selected task through the framework `design-skill` command.

## Files

- [`task_inventory.md`](/home/shuaikes/projects/agent/SciSkillUniverse/experiments/sc_skills/task_inventory.md)
- [`resource_groups.md`](/home/shuaikes/projects/agent/SciSkillUniverse/experiments/sc_skills/resource_groups.md)
- [`batch_design_manifest.json`](/home/shuaikes/projects/agent/SciSkillUniverse/experiments/sc_skills/batch_design_manifest.json)

## Batch usage

```bash
python3 scripts/batch_design_experiment_skills.py --print-only --limit 3
python3 scripts/batch_design_experiment_skills.py --task annotation
python3 scripts/batch_design_experiment_skills.py --all
python3 scripts/batch_design_experiment_skills.py --all --run-framework --verification-mode validate
python3 scripts/validate_sc_skill_experiments.py --all --json-out scratch/reviews/sc_skills_validation.json --markdown-out scratch/reviews/sc_skills_validation.md
python3 scripts/run_ref_skill_framework_batches.py run --print-plan
```

Each task is expected to include documentation plus runnable toy scaffolding under `scripts/`, `examples/`, `tests/`, and `assets/`.

