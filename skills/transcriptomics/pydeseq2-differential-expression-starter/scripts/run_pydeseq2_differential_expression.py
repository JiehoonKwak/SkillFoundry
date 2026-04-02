#!/usr/bin/env python3
"""Run a deterministic toy PyDESeq2 differential-expression workflow."""

from __future__ import annotations

import argparse
import csv
import io
import json
import warnings
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path


warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", message="invalid value encountered in slogdet")


DEFAULT_COUNTS = Path(__file__).resolve().parents[1] / "examples" / "toy_counts.tsv"
DEFAULT_METADATA = Path(__file__).resolve().parents[1] / "examples" / "toy_metadata.tsv"


def load_stack():
    try:
        import numpy as np
        import pandas as pd
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError as exc:  # pragma: no cover - exercised via runtime tests
        raise SystemExit(
            "This script requires the transcriptomics prefix at slurm/envs/transcriptomics."
        ) from exc
    return np, pd, DeseqDataSet, DeseqStats


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_inputs(counts_path: Path, metadata_path: Path):
    _, pd, _, _ = load_stack()
    count_rows = read_tsv(counts_path)
    metadata_rows = read_tsv(metadata_path)
    if not count_rows:
        raise SystemExit(f"No count rows found in {counts_path}")
    if not metadata_rows:
        raise SystemExit(f"No metadata rows found in {metadata_path}")

    genes = [key for key in count_rows[0].keys() if key != "sample"]
    counts = pd.DataFrame(
        [
            {gene: int(row[gene]) for gene in genes}
            for row in count_rows
        ],
        index=[row["sample"] for row in count_rows],
    )
    metadata = pd.DataFrame(
        [{"condition": row["condition"]} for row in metadata_rows],
        index=[row["sample"] for row in metadata_rows],
    )
    missing_samples = sorted(set(counts.index) - set(metadata.index))
    if missing_samples:
        raise SystemExit(f"Metadata is missing samples: {', '.join(missing_samples)}")
    metadata = metadata.loc[counts.index]
    return counts, metadata


def run_pydeseq2(counts_path: Path, metadata_path: Path, control: str, case: str) -> dict:
    np, pd, DeseqDataSet, DeseqStats = load_stack()
    counts, metadata = load_inputs(counts_path=counts_path, metadata_path=metadata_path)

    np.seterr(all="ignore")
    with io.StringIO() as sink, redirect_stdout(sink), redirect_stderr(sink):
        dds = DeseqDataSet(counts=counts, metadata=metadata, design="~condition", refit_cooks=True)
        dds.deseq2()
        stats = DeseqStats(dds, contrast=("condition", case, control))
        stats.summary()
    results = stats.results_df.sort_values("padj").reset_index(names="gene")

    top_results = []
    for _, row in results.head(4).iterrows():
        top_results.append(
            {
                "gene": row["gene"],
                "log2_fold_change": round(float(row["log2FoldChange"]), 6),
                "adjusted_p_value": float(row["padj"]),
            }
        )

    significant = int((results["padj"] < 0.05).sum())
    return {
        "sample_count": int(counts.shape[0]),
        "gene_count": int(counts.shape[1]),
        "design": "~condition",
        "contrast": {"factor": "condition", "case": case, "control": control},
        "significant_gene_count": significant,
        "top_gene": top_results[0]["gene"],
        "top_adjusted_p_value": top_results[0]["adjusted_p_value"],
        "top_results": top_results,
        "mean_count": round(float(pd.DataFrame(counts).to_numpy().mean()), 6),
    }


def write_json(payload: dict, out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", type=Path, default=DEFAULT_COUNTS, help="Tab-separated count matrix.")
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA, help="Tab-separated sample metadata.")
    parser.add_argument("--control", default="A", help="Control level for the contrast.")
    parser.add_argument("--case", default="B", help="Case level for the contrast.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    payload = run_pydeseq2(
        counts_path=args.counts,
        metadata_path=args.metadata,
        control=args.control,
        case=args.case,
    )
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
