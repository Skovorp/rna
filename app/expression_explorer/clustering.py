"""Sample-level dimensionality reduction for expression matrices."""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler

from expression_explorer.data import ExpressionDataset


METHODS = ("PCA", "UMAP", "t-SNE")


def sample_embedding(
    dataset: ExpressionDataset,
    method: str,
    variable_genes: int = 2_000,
) -> tuple[pd.DataFrame, str, str, str]:
    """Embed biological samples using the most variable log-TPM genes."""
    if method not in METHODS:
        raise ValueError(f"Unknown clustering method: {method}")
    if len(dataset.sample_columns) < 3:
        raise ValueError("At least three samples are required for clustering.")

    sample_by_gene = np.log2(dataset.values.to_numpy(dtype=float).T + 1.0)
    variances = sample_by_gene.var(axis=0)
    variable = np.flatnonzero(variances > 0)
    if len(variable) < 2:
        raise ValueError("At least two variable genes are required for clustering.")
    keep = variable[np.argsort(variances[variable], kind="stable")[-variable_genes:]]
    scaled = StandardScaler().fit_transform(sample_by_gene[:, keep])

    if method == "PCA":
        model = PCA(n_components=2, random_state=42)
        coordinates = model.fit_transform(scaled)
        explained = model.explained_variance_ratio_ * 100
        x_label = f"PC1 ({explained[0]:.1f}%)"
        y_label = f"PC2 ({explained[1]:.1f}%)"
        details = f"{len(keep):,} most-variable genes · PC1 + PC2 explain {explained.sum():.1f}%"
    elif method == "UMAP":
        from umap import UMAP

        neighbors = min(15, len(dataset.sample_columns) - 1)
        coordinates = UMAP(
            n_components=2,
            n_neighbors=neighbors,
            min_dist=0.2,
            metric="euclidean",
            random_state=42,
            n_jobs=1,
        ).fit_transform(scaled)
        x_label = "UMAP 1"
        y_label = "UMAP 2"
        details = f"{len(keep):,} most-variable genes · {neighbors} neighbors"
    else:
        perplexity = min(30.0, max(2.0, (len(dataset.sample_columns) - 1) / 3.0))
        coordinates = TSNE(
            n_components=2,
            perplexity=perplexity,
            init="pca",
            learning_rate="auto",
            random_state=42,
        ).fit_transform(scaled)
        x_label = "t-SNE 1"
        y_label = "t-SNE 2"
        details = f"{len(keep):,} most-variable genes · perplexity {perplexity:g}"

    metadata = (
        dataset.samples.set_index("sample", drop=False)
        .reindex(dataset.sample_columns)
        .reset_index(drop=True)
        .copy()
    )
    metadata["x"] = coordinates[:, 0]
    metadata["y"] = coordinates[:, 1]
    return metadata, x_label, y_label, details
