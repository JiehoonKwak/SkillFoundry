#!/usr/bin/env python3
"""Generate lightweight site data products from the registry."""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "registry"
SITE = ROOT / "site"
FRAMEWORK_REPORTS = ROOT / "reports" / "framework-runs"
FRAMEWORK_STATE = ROOT / "scratch" / "framework"

TAXONOMY_DOMAIN_MAP = {
    "scientific_knowledge_access_and_method_extraction": ["scientific-knowledge"],
    "data_acquisition_and_dataset_handling": ["data-acquisition-and-dataset-handling"],
    "genomics": ["genomics"],
    "transcriptomics": ["transcriptomics"],
    "epigenomics_and_chromatin": ["epigenomics-and-chromatin"],
    "proteomics_and_protein_biology": ["proteomics"],
    "metabolomics_and_other_omics": ["metabolomics-and-other-omics"],
    "structural_biology_and_molecular_modeling": ["structural-biology"],
    "systems_biology_and_network_science": ["systems-biology"],
    "imaging_and_phenotype_analysis": ["imaging-and-phenotype-analysis"],
    "drug_discovery_and_cheminformatics": ["drug-discovery-and-cheminformatics"],
    "computational_chemistry_and_molecular_simulation": ["computational-chemistry-and-molecular-simulation"],
    "materials_science_and_engineering": ["materials-science-and-engineering"],
    "earth_climate_and_geospatial_science": ["earth-climate-and-geospatial-science"],
    "clinical_biomedical_data_science": ["clinical-biomedical-data-science"],
    "statistical_and_machine_learning_foundations_for_science": ["statistical-and-machine-learning-foundations-for-science"],
    "scientific_agents_and_automation": ["scientific-agents-and-automation"],
    "reproducible_workflows_and_workflow_engines": ["reproducible-workflows"],
    "hpc_slurm_and_scaling": ["hpc"],
    "visualization_and_reporting": ["visualization-and-reporting"],
    "meta_maintenance": ["meta-maintenance"],
    "neuroscience_and_neuroimaging": ["neuroscience-and-neuroimaging"],
    "physics_and_astronomy": ["physics-and-astronomy"],
    "ecology_evolution_and_biodiversity": ["ecology-evolution-and-biodiversity"],
    "agriculture_food_and_plant_science": ["agriculture-food-and-plant-science"],
    "robotics_lab_automation_and_scientific_instrumentation": ["robotics-lab-automation-and-scientific-instrumentation"],
    "scientific_computing_and_numerical_methods": ["scientific-computing-and-numerical-methods"],
}

TAXONOMY_DISPLAY_NAMES = {
    "scientific_knowledge_access_and_method_extraction": "Scientific Knowledge Access and Method Extraction",
    "data_acquisition_and_dataset_handling": "Data Acquisition and Dataset Handling",
    "genomics": "Genomics",
    "transcriptomics": "Transcriptomics",
    "epigenomics_and_chromatin": "Epigenomics and Chromatin",
    "proteomics_and_protein_biology": "Proteomics and Protein Biology",
    "metabolomics_and_other_omics": "Metabolomics and Other Omics",
    "structural_biology_and_molecular_modeling": "Structural Biology and Molecular Modeling",
    "systems_biology_and_network_science": "Systems Biology and Network Science",
    "imaging_and_phenotype_analysis": "Imaging and Phenotype Analysis",
    "drug_discovery_and_cheminformatics": "Drug Discovery and Cheminformatics",
    "computational_chemistry_and_molecular_simulation": "Computational Chemistry and Molecular Simulation",
    "materials_science_and_engineering": "Materials Science and Engineering",
    "earth_climate_and_geospatial_science": "Earth, Climate, and Geospatial Science",
    "clinical_biomedical_data_science": "Clinical / Biomedical Data Science",
    "statistical_and_machine_learning_foundations_for_science": "Statistical and Machine Learning Foundations for Science",
    "scientific_agents_and_automation": "Scientific Agents and Automation",
    "reproducible_workflows_and_workflow_engines": "Reproducible Workflows and Workflow Engines",
    "hpc_slurm_and_scaling": "HPC, Slurm, and Scaling",
    "visualization_and_reporting": "Visualization and Reporting",
    "meta_maintenance": "Meta-maintenance",
    "neuroscience_and_neuroimaging": "Neuroscience and Neuroimaging",
    "physics_and_astronomy": "Physics and Astronomy",
    "ecology_evolution_and_biodiversity": "Ecology, Evolution, and Biodiversity",
    "agriculture_food_and_plant_science": "Agriculture, Food, and Plant Science",
    "robotics_lab_automation_and_scientific_instrumentation": "Robotics, Lab Automation, and Scientific Instrumentation",
    "scientific_computing_and_numerical_methods": "Scientific Computing and Numerical Methods",
}

TOPIC_PATH_LEAF_ALIASES = {
    ("scientific-knowledge", "method-extraction", None): "method-section-extraction",
    ("data-acquisition-and-dataset-handling", "metadata-discovery", None): "public-dataset-discovery",
    ("genomics", "gene-annotation", None): "annotation-and-effect-prediction",
    ("transcriptomics", "single-cell-rna-seq", "preprocessing"): "single-cell-rna-seq-preprocessing",
    ("transcriptomics", "single-cell-rna-seq", "integration"): "single-cell-integration-batch-correction",
    ("transcriptomics", "single-cell-rna-seq", "data-model"): "single-cell-rna-seq-preprocessing",
    ("transcriptomics", "single-cell-rna-seq", "dataset-discovery"): "multi-sample-atlas-workflows",
    ("proteomics", "dataset-discovery", None): "proteomics-dataset-discovery",
    ("proteomics", "protein-annotation", None): "protein-accession-and-metadata-lookup",
    ("structural-biology", "pdb-utilities", None): "pdb-mmcif-utilities",
    ("systems-biology", "gene-set-enrichment", "reactome"): "reactome-identifier-enrichment",
    ("systems-biology", "gene-set-enrichment", "bioconductor"): "gene-set-tooling-from-bioconductor",
    ("systems-biology", "gene-set-enrichment", "ontology"): "gene-set-enrichment",
    ("systems-biology", "pathway-knowledgebases", None): "reactome-event-lookup",
    ("drug-discovery-and-cheminformatics", "compound-search", None): "compound-identifier-lookup",
    ("drug-discovery-and-cheminformatics", "molecular-featurization", None): "fingerprints-and-similarity-search",
    ("computational-chemistry-and-molecular-simulation", "molecular-dynamics", None): "molecular-dynamics-execution",
    ("computational-chemistry-and-molecular-simulation", "quantum-chemistry", None): "single-point-energy-calculations",
    ("materials-science-and-engineering", "composition-based-featurization", None): "composition-featurization",
    ("earth-climate-and-geospatial-science", "geospatial-raster-vector-preprocessing", "geopandas"): "geospatial-feature-engineering",
    ("earth-climate-and-geospatial-science", "geospatial-raster-vector-preprocessing", "xarray"): "time-series-environmental-data-handling",
    ("clinical-biomedical-data-science", "clinical-trials", None): "cohort-extraction",
    ("reproducible-workflows", "workflow-engines", "snakemake"): "snakemake",
    ("reproducible-workflows", "workflow-engines", "nextflow"): "nextflow-nf-core",
    ("reproducible-workflows", "workflow-engines", "nf-core"): "nextflow-nf-core",
    ("reproducible-workflows", "workflow-engines", "cwl"): "cwl-command-line-tools",
    ("reproducible-workflows", "workflow-engines", "wdl"): "wdl-tasks",
    ("hpc", "slurm", "accounting"): "monitoring-and-accounting",
    ("hpc", "slurm", "batch-jobs"): "batch-jobs",
    ("hpc", "slurm", None): "batch-jobs",
}


def load_jsonl(path: Path) -> list[dict]:
    records: list[dict] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        records.append(json.loads(line))
    return records


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def read_json_if_exists(path: Path) -> dict | None:
    if not path.exists():
        return None
    return load_json(path)


def slugify_label(value: str) -> str:
    value = value.strip().lower().replace("/", " ").replace("_", " ")
    return "-".join(part for part in value.split() if part)


def resolve_leaf_slug(topic_path: list[str], valid_leaf_slugs: set[str]) -> str | None:
    if len(topic_path) < 2:
        return None
    domain = topic_path[0]
    second = topic_path[1]
    third = topic_path[2] if len(topic_path) > 2 else None
    for key in ((domain, second, third), (domain, second, None)):
        leaf_slug = TOPIC_PATH_LEAF_ALIASES.get(key)
        if leaf_slug in valid_leaf_slugs:
            return leaf_slug
    if second in valid_leaf_slugs:
        return second
    if third and third in valid_leaf_slugs:
        return third
    return None


def build_tree(skills: list[dict], resources: list[dict], taxonomy: dict) -> dict:
    root: dict = {
        "name": "SciSkillUniverse",
        "children": [],
        "skill_count": len(skills),
        "resource_count": len(resources),
        "taxonomy_domain_count": len(taxonomy),
        "covered_leaf_count": 0,
        "frontier_leaf_count": 0,
        "todo_leaf_count": 0,
    }
    skills_by_domain: dict[str, list[dict]] = {}
    resources_by_domain: dict[str, list[dict]] = {}

    for skill in skills:
        skills_by_domain.setdefault(skill["domain"], []).append(skill)

    for resource in resources:
        resources_by_domain.setdefault(resource["topic_path"][0], []).append(resource)

    for taxonomy_key, leaf_topics in taxonomy.items():
        registry_domains = TAXONOMY_DOMAIN_MAP.get(taxonomy_key, [taxonomy_key.replace("_", "-")])
        domain_skills: list[dict] = []
        domain_resources: list[dict] = []
        for domain in registry_domains:
            domain_skills.extend(skills_by_domain.get(domain, []))
            domain_resources.extend(resources_by_domain.get(domain, []))

        coverage_status = "covered" if domain_skills else ("frontier" if domain_resources else "empty")
        valid_leaf_slugs = {slugify_label(topic) for topic in leaf_topics}
        resolved_skill_leaf = {
            skill["slug"]: resolve_leaf_slug(skill.get("topic_path", []), valid_leaf_slugs)
            for skill in domain_skills
        }
        resolved_resource_leaf = {
            resource["resource_id"]: resolve_leaf_slug(resource.get("topic_path", []), valid_leaf_slugs)
            for resource in domain_resources
        }
        leaf_children = []
        for leaf_topic in leaf_topics:
            leaf_slug = slugify_label(leaf_topic)
            leaf_skills = [
                skill for skill in domain_skills if resolved_skill_leaf.get(skill["slug"]) == leaf_slug
            ]
            leaf_resources = [
                resource
                for resource in domain_resources
                if resolved_resource_leaf.get(resource["resource_id"]) == leaf_slug
            ]
            leaf_children.append(
                {
                    "name": leaf_topic,
                    "topic_slug": leaf_slug,
                    "children": [],
                    "skills": [
                    {
                        "name": skill["name"],
                        "slug": skill["slug"],
                        "status": skill["status"],
                        "path": skill["path"],
                    }
                    for skill in sorted(leaf_skills, key=lambda item: item["name"])
                ],
                "resources": [
                    {
                        "resource_id": resource["resource_id"],
                        "canonical_name": resource["canonical_name"],
                        "url": resource["url"],
                    }
                    for resource in sorted(leaf_resources, key=lambda item: item["canonical_name"])
                ],
                    "resource_count": len(leaf_resources),
                    "skill_count": len(leaf_skills),
                    "coverage_status": "covered" if leaf_skills else ("frontier" if leaf_resources else "todo"),
                    "leaf": True,
                }
            )
        covered_leaf_count = sum(1 for child in leaf_children if child["coverage_status"] == "covered")
        frontier_leaf_count = sum(1 for child in leaf_children if child["coverage_status"] == "frontier")
        todo_leaf_count = sum(1 for child in leaf_children if child["coverage_status"] == "todo")
        root["covered_leaf_count"] += covered_leaf_count
        root["frontier_leaf_count"] += frontier_leaf_count
        root["todo_leaf_count"] += todo_leaf_count
        root["children"].append(
            {
                "name": TAXONOMY_DISPLAY_NAMES.get(taxonomy_key, taxonomy_key.replace("_", " ").title()),
                "taxonomy_key": taxonomy_key,
                "registry_domains": registry_domains,
                "skill_count": len(domain_skills),
                "resource_count": len(domain_resources),
                "coverage_status": coverage_status,
                "skills": [
                    {
                        "name": skill["name"],
                        "slug": skill["slug"],
                        "status": skill["status"],
                        "path": skill["path"],
                        "topic_path": skill["topic_path"],
                    }
                    for skill in sorted(domain_skills, key=lambda item: item["name"])
                ],
                "resources": [
                    {
                        "resource_id": resource["resource_id"],
                        "canonical_name": resource["canonical_name"],
                        "status": resource["status"],
                        "url": resource["url"],
                    }
                    for resource in sorted(domain_resources, key=lambda item: item["canonical_name"])
                ],
                "children": leaf_children,
                "covered_leaf_count": covered_leaf_count,
                "frontier_leaf_count": frontier_leaf_count,
                "todo_leaf_count": todo_leaf_count,
            }
        )
    return root


def build_graph(skills: list[dict], resources: list[dict]) -> dict:
    nodes = []
    edges = []

    for resource in resources:
        nodes.append(
            {
                "id": resource["resource_id"],
                "label": resource["canonical_name"],
                "type": "resource",
                "status": resource["status"],
                "topic_path": resource["topic_path"],
                "url": resource["url"],
            }
        )

    for skill in skills:
        nodes.append(
            {
                "id": skill["slug"],
                "label": skill["name"],
                "type": "skill",
                "status": skill["status"],
                "topic_path": skill["topic_path"],
                "path": skill["path"],
            }
        )

        for resource_id in skill["source_resource_ids"]:
            edges.append({"source": skill["slug"], "target": resource_id, "type": "derived_from"})

        for related in skill.get("related_skills", []):
            edges.append({"source": skill["slug"], "target": related, "type": "related_skill"})

    return {"nodes": nodes, "edges": edges}


def parse_framework_run_timestamp(run_id: str) -> tuple[str | None, str]:
    try:
        dt = datetime.strptime("-".join(run_id.split("-")[:2]), "%Y%m%d-%H%M%S")
    except ValueError:
        return None, run_id
    return dt.isoformat(), dt.strftime("%Y-%m-%d %H:%M:%S")


def summarize_framework_stage(stage_record: dict) -> dict:
    codex = stage_record.get("codex", {})
    parsed = codex.get("parsed_message", {})
    blockers = parsed.get("blockers", []) or []
    health = "clean"
    if codex.get("returncode", 0) != 0 or blockers or parsed.get("stage") == "unparsed":
        health = "attention"
    return {
        "stage": stage_record.get("stage", "unknown"),
        "summary": parsed.get("summary", "No structured stage summary."),
        "selected_topics": parsed.get("selected_topics", []),
        "skills_touched_count": len(parsed.get("skills_touched", [])),
        "resources_touched_count": len(parsed.get("resources_touched", [])),
        "tests_run_count": len(parsed.get("tests_run", [])),
        "blocker_count": len(blockers),
        "health": health,
    }


def summarize_framework_verification(records: list[dict]) -> dict:
    passed = sum(1 for record in records if record.get("returncode") == 0)
    failed = sum(1 for record in records if record.get("returncode") != 0)
    return {
        "count": len(records),
        "passed": passed,
        "failed": failed,
        "labels": [record.get("label", "unknown") for record in records],
    }


def build_framework_runs(
    skill_count: int,
    resource_count: int,
    covered_leaf_count: int,
    frontier_leaf_count: int,
    todo_leaf_count: int,
) -> dict:
    manifests = sorted(FRAMEWORK_REPORTS.glob("*/manifest.json"), key=lambda path: path.parent.name, reverse=True)
    runs: list[dict] = []
    for manifest_path in manifests:
        manifest = load_json(manifest_path)
        run_dir = manifest_path.parent
        run_id = run_dir.name
        started_at, started_label = parse_framework_run_timestamp(run_id)
        stage_records = manifest.get("stage_results") or (
            [manifest["stage_result"]] if "stage_result" in manifest else []
        )
        stages = [summarize_framework_stage(stage) for stage in stage_records]
        verification = summarize_framework_verification(manifest.get("verification", []))
        stage_attention_count = sum(stage["health"] != "clean" for stage in stages)
        blocker_count = sum(stage["blocker_count"] for stage in stages) + verification["failed"]
        health = "clean" if blocker_count == 0 and stage_attention_count == 0 else "attention"
        latest_summary = stages[-1]["summary"] if stages else "No stage results recorded."
        final_status = manifest.get("final_status") or manifest.get("initial_status") or {}
        if isinstance(final_status, dict) and "summary" in final_status:
            final_status = final_status["summary"]
        runs.append(
            {
                "id": run_id,
                "mode": manifest.get("mode", "unknown"),
                "run_label": manifest.get("run_label", run_id),
                "run_dir": str(run_dir.relative_to(ROOT)),
                "started_at": started_at,
                "started_label": started_label,
                "loops": manifest.get("loops", 1),
                "verification_mode": manifest.get("verification_mode", "none"),
                "task_prompt": manifest.get("task_prompt"),
                "summary": latest_summary,
                "health": health,
                "blocker_count": blocker_count,
                "stage_count": len(stages),
                "stages": stages,
                "verification": verification,
                "final_status": final_status,
                "links": {
                    "manifest": str(manifest_path.relative_to(ROOT)),
                },
            }
        )
    latest_status = read_json_if_exists(FRAMEWORK_STATE / "latest_status.json")
    if isinstance(latest_status, dict):
        summary = latest_status.setdefault("summary", {})
        summary["skill_count"] = skill_count
        summary["resource_count"] = resource_count
        summary["covered_leaf_count"] = covered_leaf_count
        summary["frontier_leaf_count"] = frontier_leaf_count
        summary["todo_leaf_count"] = todo_leaf_count
    return {
        "run_count": len(runs),
        "clean_run_count": sum(1 for run in runs if run["health"] == "clean"),
        "attention_run_count": sum(1 for run in runs if run["health"] != "clean"),
        "latest_status": latest_status,
        "runs": runs,
    }


def main() -> int:
    skills = load_jsonl(REGISTRY / "skills.jsonl")
    resources = load_jsonl(REGISTRY / "resources_dedup.jsonl")
    taxonomy = load_json(REGISTRY / "taxonomy.yaml")

    skill_entries = []
    for skill in skills:
        entry = dict(skill)
        skill_path = ROOT / skill["path"]
        entry["source_resource_count"] = len(skill.get("source_resource_ids", []))
        entry["links"] = {
            "skill": f"{skill['path']}/SKILL.md",
            "metadata": f"{skill['path']}/metadata.yaml",
            "refs": f"{skill['path']}/refs.md",
            "examples": f"{skill['path']}/examples",
            "tests": f"{skill['path']}/tests",
        }
        entry["file_count"] = sum(1 for _ in skill_path.rglob("*") if _.is_file())
        skill_entries.append(entry)

    tree = build_tree(skills, resources, taxonomy)
    write_json(SITE / "skills.json", {"skills": skill_entries, "count": len(skill_entries)})
    write_json(SITE / "tree.json", tree)
    write_json(SITE / "graph.json", build_graph(skills, resources))
    framework_runs = build_framework_runs(
        skill_count=len(skills),
        resource_count=len(resources),
        covered_leaf_count=tree["covered_leaf_count"],
        frontier_leaf_count=tree["frontier_leaf_count"],
        todo_leaf_count=tree["todo_leaf_count"],
    )
    write_json(SITE / "framework_runs.json", framework_runs)
    print(
        f"Generated site data for {len(skills)} skills, {len(resources)} resources, "
        f"and {framework_runs['run_count']} framework runs."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
