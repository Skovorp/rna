from itertools import combinations
from pathlib import Path

import numpy as np

from expression_explorer.data import load_datasets
from expression_explorer.differential import (
    contrast_label,
    contrast_sample_counts,
    load_differential_contrasts,
    load_differential_results,
)


EXPRESSION_DIR = Path(__file__).resolve().parents[2] / "expression"


def test_manifests_cover_every_midgut_and_ovary_pair():
    datasets = load_datasets(EXPRESSION_DIR)
    contrasts = load_differential_contrasts(EXPRESSION_DIR)

    assert set(contrasts) == {"midgut", "elife"}
    assert len(contrasts["midgut"]) == 28
    assert len(contrasts["elife"]) == 55
    assert all(
        contrast.method == "DESeq2"
        for dataset_contrasts in contrasts.values()
        for contrast in dataset_contrasts
    )
    assert "ovary_paper" not in contrasts
    assert "crop" not in contrasts

    for dataset_key in ("midgut", "elife"):
        conditions = (
            datasets[dataset_key].samples["condition"].drop_duplicates().tolist()
        )
        expected_pairs = {
            frozenset(pair) for pair in combinations(conditions, 2)
        }
        actual_pairs = {
            frozenset((contrast.reference, contrast.target))
            for contrast in contrasts[dataset_key]
        }
        assert actual_pairs == expected_pairs
        assert all(
            contrast_sample_counts(datasets[dataset_key], contrast) == (3, 3)
            for contrast in contrasts[dataset_key]
        )


def test_differential_results_are_loaded_without_recomputing_statistics():
    dataset = load_datasets(EXPRESSION_DIR)["midgut"]
    contrast = load_differential_contrasts(EXPRESSION_DIR)["midgut"][0]
    results = load_differential_results(dataset, contrast)

    assert len(results) > 10_000
    assert results["stable_id"].is_unique
    assert results["gene"].notna().all()
    assert results["fdr"].dropna().between(0, 1).all()
    assert np.isclose(
        results["fold_change"], 2 ** results["log2_fold_change"]
    ).all()
