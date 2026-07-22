import numpy as np
import pandas as pd

from expression_explorer.comparison import benjamini_hochberg, compare_conditions
from expression_explorer.data import ExpressionDataset


def test_benjamini_hochberg_preserves_order_and_monotonicity():
    adjusted = benjamini_hochberg(np.array([0.04, 0.001, 0.03, np.nan]))

    np.testing.assert_allclose(adjusted, [0.0533333333, 0.004, 0.0533333333, 1.0])


def test_compare_conditions_reports_tpm_difference_pvalue_and_fdr():
    sample_names = [f"A{i}" for i in range(4)] + [f"B{i}" for i in range(4)]
    values = pd.DataFrame(
        [
            [1.0, 1.2, 0.9, 1.1, 20.0, 21.0, 19.0, 22.0],
            [30.0, 28.0, 31.0, 29.0, 2.0, 2.2, 1.8, 2.1],
            [5.0, 5.1, 4.9, 5.0, 5.0, 5.1, 4.9, 5.0],
        ],
        columns=sample_names,
    )
    genes = pd.DataFrame(
        {
            "row_id": [0, 1, 2],
            "display_name": ["HigherB", "HigherA", "Same"],
            "stable_id": ["ID1", "ID2", "ID3"],
        }
    )
    samples = pd.DataFrame(
        {
            "sample": sample_names,
            "group": ["A"] * 4 + ["B"] * 4,
        }
    ).set_index("sample", drop=False)
    dataset = ExpressionDataset(
        key="test",
        label="Test",
        paper="Test",
        annotation_version="Test",
        genes=genes,
        values=values,
        samples=samples,
    )

    results, samples_a, samples_b = compare_conditions(dataset, "group", "A", "B")
    indexed = results.set_index("gene")

    assert (samples_a, samples_b) == (4, 4)
    assert indexed.loc["HigherB", "mean_tpm_b"] > indexed.loc["HigherB", "mean_tpm_a"]
    assert indexed.loc["HigherB", "log2_difference"] > 0
    assert indexed.loc["HigherA", "log2_difference"] < 0
    assert indexed.loc["Same", "p_value"] == 1.0
    assert indexed.loc["HigherB", "fdr"] < 0.05
    assert indexed.loc["HigherB", "significant"]
    assert results["fdr"].between(0, 1).all()
