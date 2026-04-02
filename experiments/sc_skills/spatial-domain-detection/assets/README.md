# Assets

This directory belongs to the experiment skill `spatial-domain-detection`.

The current starter uses deterministic toy data only. No biological claim should be drawn from these assets.

## What the starter actually computes

- A SpaGCN-like spatial graph that combines coordinate distance with a lightweight histology RGB proxy.
- A GraphST-like coordinate graph, graph-diffused embedding, and deterministic k-means partition.
- Domain-vs-rest marker effect sizes on normalized toy expression.
- Machine-readable QC for graph weights, per-spot method scores, and method comparison.

## What is approximated

- Histology is represented by tiny RGB proxy vectors rather than real H&E image patches.
- SpaGCN is reduced to graph-smoothed domain signature scores instead of a trained graph convolutional network.
- GraphST is reduced to diffusion plus SVD and deterministic clustering instead of self-supervised contrastive learning.

## What remains outside starter scope

- Real-image feature extraction, large spatial matrices, GPU training, and public benchmark reproduction.
- Full Scanpy/Squidpy differential-expression or enrichment workflows on real datasets.
- Any claim that the toy partition reflects tissue biology beyond the synthetic setup.

## Upgrade path

- Replace the toy spots with a pinned public dataset slice while keeping the same deliverable names.
- Swap the RGB proxy for real histology features or pixel neighborhoods when available.
- Replace the deterministic surrogates with actual SpaGCN and GraphST package runs once the heavier dependencies are acceptable.
