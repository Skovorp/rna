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


def test_midgut_differential_manifest_has_seven_nfcore_contrasts():
    datasets = load_datasets(EXPRESSION_DIR)
    contrasts = load_differential_contrasts(EXPRESSION_DIR)

    assert set(contrasts) == {"midgut"}
    assert len(contrasts["midgut"]) == 7
    first = contrasts["midgut"][0]
    assert first.method == "DESeq2"
    assert first.reference == "female_NBF"
    assert first.target == "female_3hBF"
    assert contrast_sample_counts(datasets["midgut"], first) == (3, 3)
    assert contrast_label(datasets["midgut"], first) == (
        "Female · 3 h post-blood-meal vs Female · non-blood-fed"
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
