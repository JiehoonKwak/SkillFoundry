from __future__ import annotations

import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT = (
    ROOT
    / "skills"
    / "scientific-agents-and-automation"
    / "skill-registry-router-starter"
    / "scripts"
    / "route_skill_query.py"
)


def load_script_module():
    spec = importlib.util.spec_from_file_location("route_skill_query", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MODULE = load_script_module()


class SkillRegistryRouterTests(unittest.TestCase):
    def test_route_query_prefers_pymc_for_bayesian_regression(self) -> None:
        payload = MODULE.route_query("Bayesian regression with uncertainty intervals", top_k=3)
        self.assertGreaterEqual(payload["match_count"], 1)
        self.assertEqual(payload["matches"][0]["slug"], "pymc-bayesian-linear-regression-starter")
        matched_slugs = {match["slug"] for match in payload["matches"]}
        self.assertNotIn("matplotlib-publication-plot-starter", matched_slugs)
        self.assertNotIn("skill-registry-router-starter", matched_slugs)

    def test_cli_runtime_for_geospatial_query(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "router.json"
            subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    "--query",
                    "vector spatial join for points and polygons",
                    "--top-k",
                    "3",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=60,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["matches"][0]["slug"], "geopandas-spatial-join-starter")


if __name__ == "__main__":
    unittest.main()
