import numpy as np
from sklearn.decomposition import PCA

from analysis.compare_elife_tpm import calculate_pca


def test_comparison_pca_matches_plotpca_selection_without_gene_scaling():
    published_log = np.array(
        [
            [0.0, 0.0, 10.0, 10.0],
            [0.0, 1.0, 0.0, 1.0],
            [0.0, 2.0, 4.0, 6.0],
            [3.0, 3.0, 3.0, 3.0],
        ]
    )
    reanalysis_log = np.array(
        [
            [0.0, 1.0, 0.0, 1.0],
            [0.0, 0.0, 12.0, 12.0],
            [0.0, 3.0, 6.0, 9.0],
            [3.0, 3.0, 3.0, 3.0],
        ]
    )
    published_tpm = np.exp2(published_log) - 1.0
    reanalysis_tpm = np.exp2(reanalysis_log) - 1.0

    arrays, metrics = calculate_pca(
        published_log,
        reanalysis_log,
        min_mean_tpm=0.0,
        top_genes=2,
        published_tpm=published_tpm,
        reanalysis_tpm=reanalysis_tpm,
    )

    assert set(arrays["published_keep"]) == {0, 2}
    assert set(arrays["reanalysis_keep"]) == {1, 2}
    average_within_variance = (
        published_log.var(axis=1, ddof=1)
        + reanalysis_log.var(axis=1, ddof=1)
    ) / 2
    expected_joint = set(np.argsort(average_within_variance)[-2:])
    assert set(arrays["joint_keep"]) == expected_joint
    expected = PCA(n_components=2, random_state=42).fit_transform(
        published_log[arrays["published_keep"]].T
    )
    np.testing.assert_allclose(arrays["published_scores"], expected)
    assert metrics["per_gene_standardization"] is False
    assert metrics["protocol"] == "top_500_variable_log2_tpm_unscaled"
