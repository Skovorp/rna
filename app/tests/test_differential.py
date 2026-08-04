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
    first = contrasts["midgut"][0]
    assert first.method == "DESeq2"
    assert first.reference == "female_NBF"
    assert first.target == "female_3hBF"
    assert contrast_sample_counts(datasets["midgut"], first) == (3, 3)
    assert contrast_label(datasets["midgut"], first) == (
        "Female · 3 h post-blood-meal vs Female · non-blood-fed"
    )
    ovary_first = contrasts["elife"][0]
    assert contrast_sample_counts(datasets["elife"], ovary_first) == (3, 3)
    assert contrast_label(datasets["elife"], ovary_first) == (
        "3 hours post-blood-meal vs Non-blood-fed"
    )

    for dataset_key in ("midgut", "elife"):
        conditions = (
            datasets[dataset_key].samples["condition"].drop_duplicates().tolist()
        )
        expected_pairs = set(combinations(conditions, 2))
        actual_pairs = {
            (contrast.reference, contrast.target)
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

    assert len(results) == 15_646
    assert results["stable_id"].is_unique
    assert results["gene"].notna().all()
    assert results["fdr"].dropna().between(0, 1).all()
    aael000001 = results.set_index("stable_id").loc["AAEL000001"]
    assert np.isclose(aael000001["base_mean"], 1247.83934770134)
    assert np.isclose(aael000001["log2_fold_change"], 1.59030613602752)
    assert np.isclose(
        aael000001["fold_change"], 2 ** aael000001["log2_fold_change"]
    )
