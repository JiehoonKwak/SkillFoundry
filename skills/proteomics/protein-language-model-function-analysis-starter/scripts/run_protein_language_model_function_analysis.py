#!/usr/bin/env python3
"""Run a local protein language model embedding and function-analysis starter pipeline."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path


AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
VALID_AA = set(AA_ORDER)
HYDROPHOBIC = set("AILMFWVYC")
ACIDIC = set("DE")
BASIC = set("KRH")
POLAR = set("STNQ")
AROMATIC = set("FWY")

DEFAULT_CONFIG = {
    "backend": "mock",
    "model_family": "esm2",
    "model_id": "facebook/esm2_t33_650M_UR50D",
    "pooling": "mean",
    "embedding_dim": 32,
    "top_k_neighbors": 2,
    "label_column": "function_label",
    "sequence_id_column": "sequence_id",
}

REAL_MODEL_OPTIONS = {
    "esm2": {
        "model_id": "facebook/esm2_t33_650M_UR50D",
        "source_resource_id": "esm-transformers-docs",
        "notes": "Use raw amino-acid strings with an ESM tokenizer.",
    },
    "prott5": {
        "model_id": "Rostlab/prot_t5_xl_half_uniref50-enc",
        "source_resource_id": "prott5-model-card",
        "notes": "Space-separate amino acids and normalize uncommon residues to X.",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Protein FASTA file.")
    parser.add_argument("--labels", type=Path, default=None, help="Optional TSV/CSV with sequence_id and function_label columns.")
    parser.add_argument("--config", type=Path, default=None, help="Optional JSON config.")
    parser.add_argument("--backend", choices=["mock", "transformers"], default=None, help="Override backend from config.")
    parser.add_argument("--model-family", choices=sorted(REAL_MODEL_OPTIONS), default=None, help="High-level model family alias.")
    parser.add_argument("--model-id", default=None, help="Exact model identifier for live transformers runs.")
    parser.add_argument("--pooling", choices=["mean", "cls"], default=None, help="Embedding pooling mode.")
    parser.add_argument("--embedding-dim", type=int, default=None, help="Embedding dimension for the mock backend.")
    parser.add_argument("--embeddings-out", type=Path, required=True, help="TSV output for embeddings.")
    parser.add_argument("--summary-out", type=Path, required=True, help="JSON summary output.")
    return parser.parse_args()


def load_config(path: Path | None) -> dict:
    config = dict(DEFAULT_CONFIG)
    if path is None:
        return config
    loaded = json.loads(path.read_text(encoding="utf-8"))
    config.update(loaded)
    return config


def read_fasta(path: Path) -> list[dict[str, str]]:
    sequences: list[dict[str, str]] = []
    current_id: str | None = None
    current_lines: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id is not None:
                sequences.append({"sequence_id": current_id, "sequence": "".join(current_lines).upper()})
            current_id = line[1:].strip().split()[0]
            current_lines = []
            continue
        current_lines.append(line)
    if current_id is not None:
        sequences.append({"sequence_id": current_id, "sequence": "".join(current_lines).upper()})
    if not sequences:
        raise ValueError(f"No FASTA sequences found in {path}.")
    seen_ids: set[str] = set()
    for record in sequences:
        sequence_id = record["sequence_id"]
        sequence = record["sequence"]
        if sequence_id in seen_ids:
            raise ValueError(f"Duplicate sequence_id in FASTA: {sequence_id}")
        seen_ids.add(sequence_id)
        invalid = sorted({residue for residue in sequence if residue not in VALID_AA})
        if invalid:
            joined = ",".join(invalid)
            raise ValueError(f"Unsupported residues in {sequence_id}: {joined}")
    return sequences


def detect_delimiter(path: Path) -> str:
    header = path.read_text(encoding="utf-8").splitlines()[0]
    return "\t" if header.count("\t") >= header.count(",") else ","


def read_labels(path: Path, sequence_id_column: str, label_column: str) -> dict[str, str]:
    delimiter = detect_delimiter(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None or sequence_id_column not in reader.fieldnames or label_column not in reader.fieldnames:
            raise KeyError(
                f"Label file must include columns {sequence_id_column!r} and {label_column!r}."
            )
        labels: dict[str, str] = {}
        for row in reader:
            sequence_id = (row.get(sequence_id_column) or "").strip()
            label = (row.get(label_column) or "").strip()
            if not sequence_id or not label:
                continue
            labels[sequence_id] = label
    return labels


def normalize_vector(values: list[float]) -> list[float]:
    norm = math.sqrt(sum(value * value for value in values))
    if norm == 0.0:
        return values
    return [value / norm for value in values]


def mock_embedding(sequence: str, dim: int) -> list[float]:
    if dim < 28:
        raise ValueError("embedding_dim must be at least 28 for the mock backend.")
    length = len(sequence)
    counts = Counter(sequence)
    composition = [counts[aa] / length for aa in AA_ORDER]
    extras = [
        min(length / 128.0, 1.0),
        sum(counts[aa] for aa in HYDROPHOBIC) / length,
        sum(counts[aa] for aa in ACIDIC) / length,
        sum(counts[aa] for aa in BASIC) / length,
        sum(counts[aa] for aa in POLAR) / length,
        sum(counts[aa] for aa in AROMATIC) / length,
        (counts["G"] + counts["P"]) / length,
        (counts["C"] + counts["M"]) / length,
    ]
    hash_dims = dim - len(composition) - len(extras)
    hashed = [0.0] * hash_dims
    kmers = max(len(sequence) - 2, 1)
    for index in range(kmers):
        kmer = sequence[index:index + 3]
        digest = hashlib.sha256(kmer.encode("utf-8")).digest()
        bucket = int.from_bytes(digest[:4], "big") % hash_dims
        hashed[bucket] += 1.0 / kmers
    return normalize_vector(composition + extras + hashed)


def mean_pool(hidden_states: list[list[float]]) -> list[float]:
    if not hidden_states:
        raise ValueError("Empty hidden state list.")
    width = len(hidden_states[0])
    pooled = [0.0] * width
    for row in hidden_states:
        for index, value in enumerate(row):
            pooled[index] += value
    return [value / len(hidden_states) for value in pooled]


def transformers_embeddings(sequences: list[dict[str, str]], model_family: str, model_id: str, pooling: str) -> list[list[float]]:
    try:
        import torch
        from transformers import AutoModel, AutoTokenizer
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("The transformers backend requires both torch and transformers.") from exc

    tokenizer = AutoTokenizer.from_pretrained(model_id)
    model = AutoModel.from_pretrained(model_id)
    model.eval()
    embeddings: list[list[float]] = []
    with torch.no_grad():
        for record in sequences:
            sequence = record["sequence"]
            if model_family == "prott5":
                sequence = " ".join(sequence)
            tokens = tokenizer(sequence, return_tensors="pt")
            outputs = model(**tokens)
            hidden = outputs.last_hidden_state[0]
            if pooling == "cls":
                vector = hidden[0].tolist()
            else:
                vector = mean_pool(hidden.tolist())
            embeddings.append(normalize_vector(vector))
    return embeddings


def cosine_similarity(left: list[float], right: list[float]) -> float:
    return sum(a * b for a, b in zip(left, right))


def compute_nearest_neighbors(sequence_ids: list[str], embeddings: list[list[float]], top_k: int) -> list[dict[str, object]]:
    neighbors: list[dict[str, object]] = []
    for index, sequence_id in enumerate(sequence_ids):
        scores: list[tuple[str, float]] = []
        for other_index, other_id in enumerate(sequence_ids):
            if index == other_index:
                continue
            scores.append((other_id, cosine_similarity(embeddings[index], embeddings[other_index])))
        scores.sort(key=lambda item: item[1], reverse=True)
        neighbors.append(
            {
                "sequence_id": sequence_id,
                "neighbors": [
                    {"sequence_id": other_id, "cosine_similarity": round(score, 6)}
                    for other_id, score in scores[:top_k]
                ],
            }
        )
    return neighbors


def centroid(values: list[list[float]]) -> list[float]:
    width = len(values[0])
    totals = [0.0] * width
    for row in values:
        for index, value in enumerate(row):
            totals[index] += value
    return normalize_vector([value / len(values) for value in totals])


def classify_sequences(sequence_ids: list[str], embeddings: list[list[float]], labels: dict[str, str]) -> dict[str, object]:
    labeled_ids = [sequence_id for sequence_id in sequence_ids if sequence_id in labels]
    if not labeled_ids:
        return {
            "label_count": 0,
            "labels_seen": [],
            "holdout_accuracy": None,
            "predictions": [],
        }

    vectors_by_label: dict[str, list[list[float]]] = {}
    for sequence_id, embedding in zip(sequence_ids, embeddings):
        label = labels.get(sequence_id)
        if label is not None:
            vectors_by_label.setdefault(label, []).append(embedding)

    predictions: list[dict[str, object]] = []
    correct = 0
    for sequence_id, embedding in zip(sequence_ids, embeddings):
        observed = labels.get(sequence_id)
        if observed is None:
            continue
        label_scores: list[tuple[str, float]] = []
        for label, vectors in vectors_by_label.items():
            candidate_vectors = list(vectors)
            if label == observed and len(candidate_vectors) > 1:
                removed = False
                retained: list[list[float]] = []
                for vector in candidate_vectors:
                    if not removed and vector == embedding:
                        removed = True
                        continue
                    retained.append(vector)
                candidate_vectors = retained or candidate_vectors
            label_scores.append((label, cosine_similarity(embedding, centroid(candidate_vectors))))
        label_scores.sort(key=lambda item: item[1], reverse=True)
        predicted = label_scores[0][0]
        if predicted == observed:
            correct += 1
        predictions.append(
            {
                "sequence_id": sequence_id,
                "observed_label": observed,
                "predicted_label": predicted,
                "confidence": round(label_scores[0][1], 6),
                "scores": [
                    {"label": label, "cosine_similarity": round(score, 6)}
                    for label, score in label_scores
                ],
            }
        )
    return {
        "label_count": len(vectors_by_label),
        "labels_seen": sorted(vectors_by_label),
        "holdout_accuracy": round(correct / len(predictions), 6) if predictions else None,
        "predictions": predictions,
    }


def sequence_statistics(sequences: list[dict[str, str]]) -> dict[str, object]:
    lengths = [len(record["sequence"]) for record in sequences]
    aa_counter = Counter("".join(record["sequence"] for record in sequences))
    top_residues = [
        {"residue": residue, "count": count}
        for residue, count in aa_counter.most_common(5)
    ]
    return {
        "min_length": min(lengths),
        "max_length": max(lengths),
        "mean_length": round(sum(lengths) / len(lengths), 3),
        "top_residues": top_residues,
    }


def write_embeddings(path: Path, sequences: list[dict[str, str]], labels: dict[str, str], embeddings: list[list[float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["sequence_id", "function_label", "sequence_length"] + [
        f"dim_{index:03d}" for index in range(len(embeddings[0]))
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for record, vector in zip(sequences, embeddings):
            row = {
                "sequence_id": record["sequence_id"],
                "function_label": labels.get(record["sequence_id"], ""),
                "sequence_length": len(record["sequence"]),
            }
            row.update({fieldnames[index + 3]: round(value, 8) for index, value in enumerate(vector)})
            writer.writerow(row)


def main() -> int:
    args = parse_args()
    config = load_config(args.config)
    if args.backend is not None:
        config["backend"] = args.backend
    if args.model_family is not None:
        config["model_family"] = args.model_family
    if args.model_id is not None:
        config["model_id"] = args.model_id
    if args.pooling is not None:
        config["pooling"] = args.pooling
    if args.embedding_dim is not None:
        config["embedding_dim"] = args.embedding_dim

    sequences = read_fasta(args.input)
    labels = {}
    if args.labels is not None:
        labels = read_labels(args.labels, config["sequence_id_column"], config["label_column"])
        missing = sorted(sequence_id for sequence_id in labels if sequence_id not in {record["sequence_id"] for record in sequences})
        if missing:
            raise ValueError(f"Labels reference unknown sequence IDs: {', '.join(missing)}")

    backend = config["backend"]
    model_family = config["model_family"]
    model_id = config["model_id"]
    pooling = config["pooling"]
    if model_family not in REAL_MODEL_OPTIONS:
        raise ValueError(f"Unsupported model_family: {model_family}")

    if backend == "mock":
        embeddings = [mock_embedding(record["sequence"], int(config["embedding_dim"])) for record in sequences]
    else:
        embeddings = transformers_embeddings(sequences, model_family=model_family, model_id=model_id, pooling=pooling)

    sequence_ids = [record["sequence_id"] for record in sequences]
    neighbor_summary = compute_nearest_neighbors(sequence_ids, embeddings, int(config["top_k_neighbors"]))
    classification = classify_sequences(sequence_ids, embeddings, labels)

    write_embeddings(args.embeddings_out, sequences, labels, embeddings)

    summary = {
        "skill_slug": "protein-language-model-function-analysis-starter",
        "backend": backend,
        "model_family": model_family,
        "model_id": model_id,
        "pooling": pooling,
        "sequence_count": len(sequences),
        "embedding_dim": len(embeddings[0]),
        "sequence_statistics": sequence_statistics(sequences),
        "label_summary": {
            "label_count": classification["label_count"],
            "labels_seen": classification["labels_seen"],
        },
        "holdout_accuracy": classification["holdout_accuracy"],
        "predictions": classification["predictions"],
        "nearest_neighbors": neighbor_summary,
        "real_model_options": REAL_MODEL_OPTIONS,
        "recommended_follow_up": [
            {
                "resource_id": "esm-transformers-docs",
                "recommendation": "Use the same CLI with --backend transformers for ESM-2 embedding extraction.",
            },
            {
                "resource_id": "prott5-model-card",
                "recommendation": "Swap model_family to prott5 when encoder-only ProtT5 embeddings are a better fit.",
            },
            {
                "resource_id": "deepfri-github",
                "recommendation": "Feed the emitted embedding or label summary into downstream sequence-to-function benchmarking or annotation workflows.",
            },
        ],
        "caveats": [
            "The bundled mock backend is for deterministic smoke testing only.",
            "Toy labels illustrate pipeline wiring and should not be interpreted as validated biology.",
        ],
    }
    args.summary_out.parent.mkdir(parents=True, exist_ok=True)
    args.summary_out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
