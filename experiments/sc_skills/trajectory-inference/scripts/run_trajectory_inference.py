#!/usr/bin/env python3
from __future__ import annotations

import argparse
import heapq
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


SKILL_DIR = Path(__file__).resolve().parents[1]
QC_ARTIFACTS = [
    "qc_preprocessing.json",
    "qc_knn_graph.tsv",
    "qc_velocity_surrogate.tsv",
    "qc_transition_matrix.tsv",
    "run_summary.json",
]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_metadata(skill_dir: Path) -> dict[str, Any]:
    return load_json(skill_dir / "metadata.yaml")


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_markdown(path: Path, title: str, sections: list[dict[str, Any]]) -> None:
    lines = [f"# {title}", ""]
    for section in sections:
        lines.append(f"## {section['name']}")
        lines.append("")
        for paragraph in section.get("paragraphs", []):
            lines.append(paragraph)
            lines.append("")
        for bullet in section.get("bullets", []):
            lines.append(f"- {bullet}")
        lines.append("")
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def validate_markdown_sections(path: Path, required_sections: list[str]) -> None:
    text = path.read_text(encoding="utf-8")
    for section in required_sections:
        heading = f"## {section}"
        if heading not in text:
            raise AssertionError(f"Missing markdown section {heading} in {path}")


def cosine_similarity(left: np.ndarray, right: np.ndarray) -> float:
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    if denominator <= 1e-12:
        return 0.0
    return float(np.dot(left, right) / denominator)


def connected_components(adjacency: np.ndarray) -> list[list[int]]:
    remaining = set(range(len(adjacency)))
    components: list[list[int]] = []
    while remaining:
        start = remaining.pop()
        stack = [start]
        component = [start]
        while stack:
            node = stack.pop()
            neighbors = np.flatnonzero(adjacency[node] > 0.0)
            for neighbor in neighbors:
                if int(neighbor) in remaining:
                    remaining.remove(int(neighbor))
                    stack.append(int(neighbor))
                    component.append(int(neighbor))
        components.append(sorted(component))
    return components


def dijkstra_distances(edge_lengths: np.ndarray, root_index: int) -> np.ndarray:
    distances = np.full(edge_lengths.shape[0], np.inf, dtype=float)
    distances[root_index] = 0.0
    heap: list[tuple[float, int]] = [(0.0, root_index)]

    while heap:
        current_distance, current = heapq.heappop(heap)
        if current_distance > distances[current]:
            continue
        for neighbor in np.flatnonzero(np.isfinite(edge_lengths[current])):
            step = float(edge_lengths[current, neighbor])
            candidate = current_distance + step
            if candidate + 1e-12 < distances[neighbor]:
                distances[neighbor] = candidate
                heapq.heappush(heap, (candidate, int(neighbor)))
    return distances


def build_knn_graph(
    cell_ids: list[str],
    coordinates: np.ndarray,
    knn_k: int,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]], list[list[int]]]:
    if knn_k < 1:
        raise AssertionError("knn_k must be at least 1.")

    pairwise = np.linalg.norm(coordinates[:, None, :] - coordinates[None, :, :], axis=2)
    np.fill_diagonal(pairwise, np.inf)
    directed = np.zeros_like(pairwise)

    for source_index in range(len(cell_ids)):
        neighbor_indices = np.argsort(pairwise[source_index])[: min(knn_k, len(cell_ids) - 1)]
        for target_index in neighbor_indices:
            directed[source_index, target_index] = 1.0 / float(pairwise[source_index, target_index])

    adjacency = np.maximum(directed, directed.T)
    edge_lengths = np.full_like(pairwise, np.inf)
    mask = adjacency > 0.0
    edge_lengths[mask] = pairwise[mask]
    np.fill_diagonal(edge_lengths, 0.0)

    rows: list[dict[str, Any]] = []
    for source_index, source_cell in enumerate(cell_ids):
        neighbors = np.flatnonzero(adjacency[source_index] > 0.0)
        total_weight = float(adjacency[source_index, neighbors].sum())
        for target_index in neighbors:
            rows.append(
                {
                    "source_cell": source_cell,
                    "target_cell": cell_ids[int(target_index)],
                    "distance": round(float(pairwise[source_index, target_index]), 6),
                    "raw_connectivity": round(float(adjacency[source_index, target_index]), 6),
                    "row_normalized_connectivity": round(
                        float(adjacency[source_index, target_index] / total_weight),
                        6,
                    ),
                    "mutual_edge": bool(directed[source_index, target_index] > 0.0 and directed[target_index, source_index] > 0.0),
                }
            )

    return adjacency, edge_lengths, rows, connected_components(adjacency)


def compute_pseudotime(
    edge_lengths: np.ndarray,
    root_index: int,
) -> tuple[np.ndarray, np.ndarray]:
    root_distances = dijkstra_distances(edge_lengths, root_index)
    if not np.all(np.isfinite(root_distances)):
        raise AssertionError("k-NN graph must remain connected to compute pseudotime.")
    scale = float(root_distances.max())
    if scale <= 0.0:
        raise AssertionError("Pseudotime scaling failed because all cells are at zero distance.")
    pseudotime = root_distances / scale
    return root_distances, pseudotime


def build_forward_neighbors(adjacency: np.ndarray, pseudotime: np.ndarray, forward_eps: float) -> dict[int, list[int]]:
    forward_neighbors: dict[int, list[int]] = {}
    for source_index in range(len(adjacency)):
        candidates = [
            int(target_index)
            for target_index in np.flatnonzero(adjacency[source_index] > 0.0)
            if pseudotime[target_index] > pseudotime[source_index] + forward_eps
        ]
        forward_neighbors[source_index] = sorted(candidates, key=lambda index: float(pseudotime[index]))
    return forward_neighbors


def compute_velocity_surrogate(
    *,
    cell_ids: list[str],
    coordinates: np.ndarray,
    spliced: np.ndarray,
    unspliced: np.ndarray,
    adjacency: np.ndarray,
    pseudotime: np.ndarray,
    terminal_indices: set[int],
    alignment_floor: float,
    forward_eps: float,
) -> tuple[list[dict[str, Any]], np.ndarray, dict[int, list[int]]]:
    forward_neighbors = build_forward_neighbors(adjacency, pseudotime, forward_eps)
    velocity_programs = unspliced - spliced
    velocity_rows: list[dict[str, Any]] = []
    transition_scores = np.zeros_like(adjacency)

    for source_index, source_cell in enumerate(cell_ids):
        is_terminal = source_index in terminal_indices
        candidates = forward_neighbors[source_index]
        if is_terminal or not candidates:
            velocity_rows.append(
                {
                    "cell_id": source_cell,
                    "velocity_program_alpha": round(float(velocity_programs[source_index, 0]), 6),
                    "velocity_program_beta": round(float(velocity_programs[source_index, 1]), 6),
                    "velocity_dx": 0.0,
                    "velocity_dy": 0.0,
                    "velocity_magnitude": 0.0,
                    "mean_alignment": 0.0,
                    "forward_neighbor_count": 0,
                    "forward_neighbor_ids": "",
                    "is_terminal": bool(is_terminal),
                }
            )
            continue

        raw_weights: list[float] = []
        alignments: list[float] = []
        displacements: list[np.ndarray] = []
        for target_index in candidates:
            expr_delta = spliced[target_index] - spliced[source_index]
            alignment = max(cosine_similarity(velocity_programs[source_index], expr_delta), 0.0)
            pseudotime_gain = float(pseudotime[target_index] - pseudotime[source_index])
            raw_score = float(adjacency[source_index, target_index]) * max(pseudotime_gain, 0.0) * (alignment_floor + alignment)
            raw_weights.append(raw_score)
            alignments.append(alignment)
            displacements.append(coordinates[target_index] - coordinates[source_index])

        weight_array = np.asarray(raw_weights, dtype=float)
        weight_array = weight_array / weight_array.sum()
        displacement_matrix = np.vstack(displacements)
        coordinate_velocity = (weight_array[:, None] * displacement_matrix).sum(axis=0)
        transition_scores[source_index, candidates] = weight_array

        velocity_rows.append(
            {
                "cell_id": source_cell,
                "velocity_program_alpha": round(float(velocity_programs[source_index, 0]), 6),
                "velocity_program_beta": round(float(velocity_programs[source_index, 1]), 6),
                "velocity_dx": round(float(coordinate_velocity[0]), 6),
                "velocity_dy": round(float(coordinate_velocity[1]), 6),
                "velocity_magnitude": round(float(np.linalg.norm(coordinate_velocity)), 6),
                "mean_alignment": round(float(np.dot(weight_array, np.asarray(alignments, dtype=float))), 6),
                "forward_neighbor_count": len(candidates),
                "forward_neighbor_ids": ";".join(cell_ids[target_index] for target_index in candidates),
                "is_terminal": False,
            }
        )

    return velocity_rows, transition_scores, forward_neighbors


def build_transition_matrix(
    *,
    cell_ids: list[str],
    pseudotime: np.ndarray,
    transition_scores: np.ndarray,
    forward_neighbors: dict[int, list[int]],
    terminal_indices: set[int],
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    transition = np.zeros_like(transition_scores)
    rows: list[dict[str, Any]] = []

    for source_index, source_cell in enumerate(cell_ids):
        if source_index in terminal_indices:
            transition[source_index, source_index] = 1.0
            rows.append(
                {
                    "source_cell": source_cell,
                    "target_cell": source_cell,
                    "transition_probability": 1.0,
                    "pseudotime_gain": 0.0,
                    "is_absorbing": True,
                }
            )
            continue

        candidates = forward_neighbors[source_index]
        if not candidates:
            raise AssertionError(f"Cell {source_cell} has no forward neighbors in the toy graph.")

        weights = transition_scores[source_index, candidates]
        if np.any(weights <= 0.0):
            raise AssertionError(f"Cell {source_cell} produced a non-positive transition weight.")
        weights = weights / weights.sum()
        transition[source_index, candidates] = weights
        for target_index, probability in zip(candidates, weights, strict=True):
            rows.append(
                {
                    "source_cell": source_cell,
                    "target_cell": cell_ids[target_index],
                    "transition_probability": round(float(probability), 6),
                    "pseudotime_gain": round(float(pseudotime[target_index] - pseudotime[source_index]), 6),
                    "is_absorbing": False,
                }
            )

    row_sums = transition.sum(axis=1)
    if not np.allclose(row_sums, 1.0, atol=1e-6):
        raise AssertionError("Transition matrix rows must sum to 1.")

    return transition, rows


def compute_absorption_probabilities(
    transition: np.ndarray,
    terminal_indices: list[int],
) -> np.ndarray:
    all_indices = list(range(len(transition)))
    transient_indices = [index for index in all_indices if index not in terminal_indices]
    probabilities = np.zeros((len(transition), len(terminal_indices)), dtype=float)

    if transient_indices:
        q_matrix = transition[np.ix_(transient_indices, transient_indices)]
        r_matrix = transition[np.ix_(transient_indices, terminal_indices)]
        identity = np.eye(len(transient_indices), dtype=float)
        solved = np.linalg.solve(identity - q_matrix, r_matrix)
        probabilities[np.asarray(transient_indices, dtype=int), :] = solved

    for terminal_position, terminal_index in enumerate(terminal_indices):
        probabilities[terminal_index, terminal_position] = 1.0

    row_sums = probabilities.sum(axis=1)
    if not np.allclose(row_sums, 1.0, atol=1e-6):
        raise AssertionError("Absorption probabilities must sum to 1 per cell.")
    return probabilities


def entropy_from_probabilities(probabilities: np.ndarray) -> np.ndarray:
    safe = np.clip(probabilities, 1e-12, 1.0)
    entropy = -(safe * np.log2(safe)).sum(axis=1)
    return entropy


def build_outputs(payload: dict[str, Any], *, outdir: Path) -> dict[str, Any]:
    outdir.mkdir(parents=True, exist_ok=True)

    cells = payload["cells"]
    cell_ids = [cell["cell_id"] for cell in cells]
    coordinates = np.asarray([[float(cell["x"]), float(cell["y"])] for cell in cells], dtype=float)
    spliced = np.asarray([cell["spliced"] for cell in cells], dtype=float)
    unspliced = np.asarray([cell["unspliced"] for cell in cells], dtype=float)
    parameters = payload["parameters"]
    root_cell_id = payload["root_cell_id"]
    terminal_states = payload["terminal_states"]

    if len(set(cell_ids)) != len(cell_ids):
        raise AssertionError("Cell identifiers must be unique.")
    if spliced.shape != unspliced.shape:
        raise AssertionError("Spliced and unspliced arrays must have matching shapes.")
    if spliced.shape[1] != len(payload["program_names"]):
        raise AssertionError("Program names must match the length of spliced/unspliced vectors.")
    if root_cell_id not in cell_ids:
        raise AssertionError(f"Root cell {root_cell_id!r} is not present in the toy input.")

    root_index = cell_ids.index(root_cell_id)
    terminal_index_by_state: dict[str, int] = {}
    for state in terminal_states:
        state_cell_id = state["cell_id"]
        if state_cell_id not in cell_ids:
            raise AssertionError(f"Terminal cell {state_cell_id!r} is missing from the toy input.")
        terminal_index_by_state[state["state"]] = cell_ids.index(state_cell_id)

    adjacency, edge_lengths, graph_rows, components = build_knn_graph(
        cell_ids=cell_ids,
        coordinates=coordinates,
        knn_k=int(parameters["knn_k"]),
    )
    root_distances, pseudotime = compute_pseudotime(edge_lengths=edge_lengths, root_index=root_index)
    velocity_rows, transition_scores, forward_neighbors = compute_velocity_surrogate(
        cell_ids=cell_ids,
        coordinates=coordinates,
        spliced=spliced,
        unspliced=unspliced,
        adjacency=adjacency,
        pseudotime=pseudotime,
        terminal_indices=set(terminal_index_by_state.values()),
        alignment_floor=float(parameters["alignment_floor"]),
        forward_eps=float(parameters["forward_eps"]),
    )
    ordered_terminal_states = list(terminal_index_by_state)
    ordered_terminal_indices = [terminal_index_by_state[state] for state in ordered_terminal_states]
    transition_matrix, transition_rows = build_transition_matrix(
        cell_ids=cell_ids,
        pseudotime=pseudotime,
        transition_scores=transition_scores,
        forward_neighbors=forward_neighbors,
        terminal_indices=set(ordered_terminal_indices),
    )
    absorption = compute_absorption_probabilities(transition_matrix, ordered_terminal_indices)
    entropy = entropy_from_probabilities(absorption)

    lineage_confidence_threshold = float(parameters["lineage_confidence_threshold"])
    trunk_candidate_indices = np.flatnonzero(absorption.max(axis=1) < lineage_confidence_threshold)
    if len(trunk_candidate_indices) > 0:
        branchpoint_index = int(max(trunk_candidate_indices, key=lambda index: float(pseudotime[index])))
    else:
        branchpoint_index = int(np.argmin(np.abs(absorption[:, 0] - absorption[:, 1])))
    branchpoint_cell = cell_ids[branchpoint_index]

    coordinate_rows: list[dict[str, Any]] = []
    fate_rows: list[dict[str, Any]] = []
    for cell_index, cell_id in enumerate(cell_ids):
        probability_row = absorption[cell_index]
        dominant_index = int(np.argmax(probability_row))
        dominant_terminal = ordered_terminal_states[dominant_index]
        dominant_probability = float(probability_row[dominant_index])
        if cell_index == root_index:
            lineage_label = "root"
        elif dominant_probability < lineage_confidence_threshold:
            lineage_label = "trunk"
        else:
            lineage_label = dominant_terminal

        coordinate_rows.append(
            {
                "cell_id": cell_id,
                "x_coord": round(float(coordinates[cell_index, 0]), 6),
                "y_coord": round(float(coordinates[cell_index, 1]), 6),
                "root_distance": round(float(root_distances[cell_index]), 6),
                "pseudotime": round(float(pseudotime[cell_index]), 6),
                "lineage_label": lineage_label,
                "dominant_terminal_state": dominant_terminal,
                "dominant_terminal_probability": round(dominant_probability, 6),
                "entropy": round(float(entropy[cell_index]), 6),
            }
        )

        for terminal_index, terminal_state in enumerate(ordered_terminal_states):
            fate_rows.append(
                {
                    "cell_id": cell_id,
                    "terminal_state": terminal_state,
                    "probability": round(float(probability_row[terminal_index]), 6),
                }
            )

    trajectory_frame = pd.DataFrame(coordinate_rows)
    fate_frame = pd.DataFrame(fate_rows)
    graph_frame = pd.DataFrame(graph_rows)
    velocity_frame = pd.DataFrame(velocity_rows)
    transition_frame = pd.DataFrame(transition_rows)

    trajectory_frame.to_csv(outdir / "trajectory_coordinates.tsv", sep="\t", index=False)
    fate_frame.to_csv(outdir / "fate_probabilities.tsv", sep="\t", index=False)
    graph_frame.to_csv(outdir / "qc_knn_graph.tsv", sep="\t", index=False)
    velocity_frame.to_csv(outdir / "qc_velocity_surrogate.tsv", sep="\t", index=False)
    transition_frame.to_csv(outdir / "qc_transition_matrix.tsv", sep="\t", index=False)

    preprocessing = {
        "branchpoint_cell_id": branchpoint_cell,
        "cell_count": len(cell_ids),
        "component_count": len(components),
        "components": [[cell_ids[index] for index in component] for component in components],
        "connected": len(components) == 1,
        "knn_k": int(parameters["knn_k"]),
        "program_names": payload["program_names"],
        "root_cell_id": root_cell_id,
        "terminal_states": terminal_states,
    }
    write_json(outdir / "qc_preprocessing.json", preprocessing)

    lineage_counts = trajectory_frame["lineage_label"].value_counts().to_dict()
    alpha_terminal = ordered_terminal_states[0]
    beta_terminal = ordered_terminal_states[1]
    report_sections = [
        {
            "name": "Run context",
            "bullets": [
                f"Run label: {payload['run_label']}.",
                f"Computed four starter steps on {len(cell_ids)} synthetic cells: k-NN graph, shortest-path pseudotime, neighbor-delta velocity surrogate, and absorbing-chain fate probabilities.",
                f"Preprocessing checkpoints passed: unique cell ids, root cell {root_cell_id}, two absorbing terminals, and a connected graph with {len(components)} component.",
            ],
        },
        {
            "name": "Coordinates",
            "bullets": [
                f"Pseudotime spans {trajectory_frame['pseudotime'].min():.3f} to {trajectory_frame['pseudotime'].max():.3f}.",
                f"Branchpoint candidate {branchpoint_cell} sits between the two terminal states with near-balanced fate scores.",
                f"Lineage labels summarize one root cell, {lineage_counts.get('trunk', 0)} trunk cells, and branch-biased cells toward {alpha_terminal} or {beta_terminal}.",
            ],
        },
        {
            "name": "Fate summaries",
            "bullets": [
                f"{alpha_terminal} is dominant for {int((trajectory_frame['dominant_terminal_state'] == alpha_terminal).sum())} cells; {beta_terminal} is dominant for {int((trajectory_frame['dominant_terminal_state'] == beta_terminal).sum())} cells.",
                f"Terminal cells stay absorbing with probability 1.0 for their own state, while earlier trunk cells remain mixed before the branchpoint.",
                f"Velocity directions are summarized in qc_velocity_surrogate.tsv and used only as a lightweight directional cue for transition weighting.",
            ],
        },
        {
            "name": "Caveats",
            "bullets": [
                "This starter computes graph pseudotime exactly on the toy k-NN graph, but RNA velocity is only a two-program surrogate derived from synthetic spliced and unspliced scores.",
                "Fate probabilities come from a tiny absorbing Markov chain over local transitions, not from a full CellRank kernel stack or uncertainty propagation model.",
                "The outputs validate the portable experiment contract and interpretation checkpoints, not a biological claim on real single-cell data.",
            ],
        },
    ]
    write_markdown(outdir / "trajectory_report.md", "Trajectory Report", report_sections)

    summary = {
        "branchpoint_cell_id": branchpoint_cell,
        "connected_graph": bool(len(components) == 1),
        "method_steps": [
            "knn_graph_construction",
            "shortest_path_pseudotime",
            "neighbor_delta_velocity_surrogate",
            "absorbing_markov_fate_probabilities",
        ],
        "outdir": str(outdir.resolve()),
        "prediction_counts": {key: int(value) for key, value in lineage_counts.items()},
        "root_cell_id": root_cell_id,
        "terminal_states": ordered_terminal_states,
        "written_files": sorted(
            [
                "fate_probabilities.tsv",
                "qc_knn_graph.tsv",
                "qc_preprocessing.json",
                "qc_transition_matrix.tsv",
                "qc_velocity_surrogate.tsv",
                "run_summary.json",
                "trajectory_coordinates.tsv",
                "trajectory_report.md",
            ]
        ),
    }
    write_json(outdir / "run_summary.json", summary)
    return summary


def validate_outputs(skill_dir: Path, outdir: Path) -> None:
    metadata = load_metadata(skill_dir)
    for deliverable in metadata["deliverables"]:
        path = outdir / deliverable["path"]
        if not path.exists():
            raise AssertionError(f"Missing deliverable: {path}")
        if deliverable["kind"] == "tsv":
            frame = pd.read_csv(path, sep="\t")
            missing = [column for column in deliverable.get("required_columns", []) if column not in frame.columns]
            if missing:
                raise AssertionError(f"Missing TSV columns in {path}: {missing}")
            if frame.empty:
                raise AssertionError(f"Empty TSV deliverable: {path}")
        elif deliverable["kind"] == "md":
            validate_markdown_sections(path, deliverable.get("required_sections", []))
        else:
            raise AssertionError(f"Unsupported deliverable kind {deliverable['kind']!r}.")

    trajectory = pd.read_csv(outdir / "trajectory_coordinates.tsv", sep="\t")
    if trajectory["cell_id"].duplicated().any():
        raise AssertionError("Trajectory output must keep unique cell identifiers.")
    if not ((trajectory["pseudotime"] >= 0.0) & (trajectory["pseudotime"] <= 1.0)).all():
        raise AssertionError("Pseudotime values must stay within [0, 1].")

    fate = pd.read_csv(outdir / "fate_probabilities.tsv", sep="\t")
    grouped = fate.groupby("cell_id")["probability"].sum()
    if not np.allclose(grouped.to_numpy(dtype=float), 1.0, atol=1e-6):
        raise AssertionError("Fate probabilities must sum to 1 per cell.")

    for artifact in metadata.get("qc_artifacts", QC_ARTIFACTS):
        if not (outdir / artifact).exists():
            raise AssertionError(f"Missing QC artifact: {outdir / artifact}")

    preprocessing = load_json(outdir / "qc_preprocessing.json")
    if not preprocessing.get("connected", False):
        raise AssertionError("Toy graph must stay connected.")

    transition = pd.read_csv(outdir / "qc_transition_matrix.tsv", sep="\t")
    transition_sums = transition.groupby("source_cell")["transition_probability"].sum()
    if not np.allclose(transition_sums.to_numpy(dtype=float), 1.0, atol=1e-6):
        raise AssertionError("Transition matrix rows must sum to 1.")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=SKILL_DIR / "examples" / "toy_input.json")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)

    summary = build_outputs(load_json(args.input), outdir=args.outdir.resolve())
    validate_outputs(SKILL_DIR, args.outdir.resolve())
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
