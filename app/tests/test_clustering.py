import numpy as np
import pandas as pd

from expression_explorer.clustering import METHODS, sample_embedding
from expression_explorer.data import ExpressionDataset


def synthetic_dataset() -> ExpressionDataset:
    rng = np.random.default_rng(42)
    sample_names = [f"sample_{index}" for index in range(12)]
    values = rng.gamma(shape=2.0, scale=4.0, size=(120, 12))
    values[:20, 6:] += 15.0
    samples = pd.DataFrame(
        {
            "sample": sample_names,
            "condition_label": ["A"] * 6 + ["B"] * 6,
            "tissue": ["antenna"] * 12,
            "sex": ["female"] * 12,
        }
    ).set_index("sample", drop=False)
    genes = pd.DataFrame(
        {
            "row_id": np.arange(120),
            "display_name": [f"gene_{index}" for index in range(120)],
            "stable_id": [f"ID{index}" for index in range(120)],
        }
    )
    return ExpressionDataset(
        key="synthetic",
        label="Synthetic",
        paper="Test",
        annotation_version="Test",
        genes=genes,
        values=pd.DataFrame(values, columns=sample_names),
        samples=samples,
    )


def test_pca_embedding_is_sample_aligned_and_reproducible():
    dataset = synthetic_dataset()
    first, x_label, y_label, details = sample_embedding(dataset, "PCA", 80)
    second, *_ = sample_embedding(dataset, "PCA", 80)

    assert first["sample"].tolist() == dataset.sample_columns
    assert first[["x", "y"]].shape == (12, 2)
    assert np.isfinite(first[["x", "y"]]).all(axis=None)
    np.testing.assert_allclose(first[["x", "y"]], second[["x", "y"]])
    assert x_label.startswith("PC1 (")
    assert y_label.startswith("PC2 (")
    assert "80 most-variable genes" in details


def test_nonlinear_embeddings_return_two_finite_coordinates():
    dataset = synthetic_dataset()
    for method in ("UMAP", "t-SNE"):
        embedding, x_label, y_label, details = sample_embedding(dataset, method, 80)
        assert embedding[["x", "y"]].shape == (12, 2)
        assert np.isfinite(embedding[["x", "y"]]).all(axis=None)
        assert x_label.startswith(method)
        assert y_label.endswith("2")
        assert "80 most-variable genes" in details


def test_supported_cluster_methods_are_stable():
    assert METHODS == ("PCA", "UMAP", "t-SNE")
