# Suggested Checks

- Run `python3 skills/reproducible-workflows/nf-core-pipeline-list/scripts/list_nfcore_pipelines.py --limit 3`.
- Confirm that the JSON summary contains at least one `remote_workflows` entry.
- Confirm that the wrapper succeeds even though the raw nf-core CLI JSON is malformed in this environment.
