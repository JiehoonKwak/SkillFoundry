from __future__ import annotations

import json
import unittest
from pathlib import Path

import scripts.build_site as build_site


ROOT = Path(__file__).resolve().parents[2]
TAXONOMY = json.loads((ROOT / "registry" / "taxonomy.yaml").read_text(encoding="utf-8"))


class BuildSiteAliasResolutionTests(unittest.TestCase):
    def test_specific_three_part_alias_beats_generic_second_level_leaf(self) -> None:
        valid = {build_site.slugify_label(topic) for topic in TAXONOMY["systems_biology_and_network_science"]}
        self.assertEqual(
            build_site.resolve_leaf_slug(["systems-biology", "gene-set-enrichment", "reactome"], valid),
            "reactome-identifier-enrichment",
        )
        self.assertEqual(
            build_site.resolve_leaf_slug(["systems-biology", "gene-set-enrichment", "bioconductor"], valid),
            "gene-set-tooling-from-bioconductor",
        )


if __name__ == "__main__":
    unittest.main()
