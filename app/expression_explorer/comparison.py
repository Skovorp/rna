"""Exploratory condition comparisons for published TPM matrices."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

from expression_explorer.data import ExpressionDataset


def benjamini_hochberg(p_values: np.ndarray | pd.Series) -> np.ndarray:
    """Return Benjamini-Hochberg adjusted p-values in the original order."""
    values = np.asarray(p_values, dtype=float)
    values = np.where(np.isfinite(values), values, 1.0)
    values = np.clip(values, 0.0, 1.0)
    count = len(values)
    if count == 0:
        return values

    order = np.argsort(values, kind="stable")
    ranked = values[order]
    adjusted = ranked * count / np.arange(1, count + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty(count, dtype=float)
    restored[order] = np.clip(adjusted, 0.0, 1.0)
    return restored


def compare_conditions(
    dataset: ExpressionDataset,
    field: str,
    condition_a: str,
    condition_b: str,
) -> tuple[pd.DataFrame, int, int]:
    """Compare all genes between two groups with Welch tests on log2(TPM + 1).

    The returned log difference is mean(log2(TPM + 1)) in B minus A. Raw
    p-values are adjusted across all tested genes using Benjamini-Hochberg FDR.
    """
    if field not in dataset.samples.columns:
        raise ValueError(f"Unknown sample grouping field: {field}")

    metadata = dataset.samples.copy()
    metadata[field] = metadata[field].fillna("").astype(str)
    samples_a = [
        sample
        for sample in metadata.loc[metadata[field].eq(str(condition_a)), "sample"]
        if sample in dataset.values.columns
    ]
    samples_b = [
        sample
        for sample in metadata.loc[metadata[field].eq(str(condition_b)), "sample"]
        if sample in dataset.values.columns
    ]
    if len(samples_a) < 2 or len(samples_b) < 2:
        raise ValueError("Each condition needs at least two biological samples.")

    tpm_a = dataset.values[samples_a].to_numpy(dtype=float)
    tpm_b = dataset.values[samples_b].to_numpy(dtype=float)
    log_a = np.log2(tpm_a + 1.0)
    log_b = np.log2(tpm_b + 1.0)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        test = ttest_ind(log_b, log_a, axis=1, equal_var=False, nan_policy="omit")
    p_values = np.asarray(test.pvalue, dtype=float)
    p_values = np.where(np.isfinite(p_values), p_values, 1.0)

    annotations = (
        dataset.genes.drop_duplicates("row_id")
        .set_index("row_id")
        .reindex(dataset.values.index)
    )
    results = pd.DataFrame(
        {
            "gene": annotations["display_name"].fillna("").to_numpy(),
            "stable_id": annotations["stable_id"].fillna("").to_numpy(),
            "mean_tpm_a": tpm_a.mean(axis=1),
            "mean_tpm_b": tpm_b.mean(axis=1),
            "median_tpm_a": np.median(tpm_a, axis=1),
            "median_tpm_b": np.median(tpm_b, axis=1),
            "log2_difference": log_b.mean(axis=1) - log_a.mean(axis=1),
            "p_value": p_values,
            "fdr": benjamini_hochberg(p_values),
        }
    )
    results["significant"] = results["fdr"] < 0.05
    results["absolute_log2_difference"] = results["log2_difference"].abs()
    results = results.sort_values(
        ["fdr", "p_value", "absolute_log2_difference"],
        ascending=[True, True, False],
        kind="stable",
    ).reset_index(drop=True)
    return results, len(samples_a), len(samples_b)
