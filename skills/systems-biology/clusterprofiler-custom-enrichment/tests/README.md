# Tests

`test_clusterprofiler_custom_enrichment.py` always validates `--describe-toy`.
Set `BIOC_SKILL_R_LIB` to a library containing `clusterProfiler` to run the package-backed integration branch.
In this environment, `clusterProfiler` installation is blocked by missing Cairo development headers required by `gdtools`.
