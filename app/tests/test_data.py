from pathlib import Path

from expression_explorer.data import (
    classify_family,
    expression_long,
    family_members,
    load_datasets,
    load_nfcore_dataset,
    search_genes,
)


EXPRESSION_DIR = Path(__file__).resolve().parents[2] / "expression"


def test_family_classifier_accepts_numbered_gene_copies():
    for symbol in ("Ir31a1", "Ir100a.1", "AaegIr100a.2"):
        canonical = symbol.removeprefix("Aaeg")
        assert classify_family(canonical) == "Ionotropic receptors (IR)"
    assert classify_family("Ir") == "Other"


def test_dataset_dimensions_and_sample_metadata():
    datasets = load_datasets(EXPRESSION_DIR)
    assert set(datasets) >= {
        "ovary_paper",
        "elife",
        "neuro_ru",
        "neuro_legacy",
        "midgut",
        "crop",
    }
    assert datasets["neuro_ru"].values.shape == (16_176, 122)
    assert datasets["neuro_legacy"].values.shape == (17_478, 122)
    assert datasets["ovary_paper"].values.shape[1] == 33
    assert datasets["elife"].values.shape[1] == 33
    assert datasets["midgut"].values.shape[1] == 24
    assert datasets["crop"].values.shape[1] == 3
    for dataset in datasets.values():
        assert dataset.samples.index.tolist() == dataset.sample_columns
        assert dataset.samples["sample"].notna().all()

    midgut = datasets["midgut"]
    assert midgut.samples["tissue"].eq("midgut").all()
    assert midgut.samples["condition_label"].nunique() == 8
    assert not midgut.samples["condition_label"].eq("Unspecified").any()
    assert midgut.samples["condition_label"].drop_duplicates().tolist() == [
        "Non-blood-fed",
        "3 hours post-blood-meal",
        "6 hours post-blood-meal",
        "12 hours post-blood-meal",
        "24 hours post-blood-meal",
        "48 hours post-blood-meal",
        "72 hours post-blood-meal",
        "Male · non-blood-fed",
    ]
    assert midgut.samples["sex"].drop_duplicates().tolist() == ["female", "male"]

    elife = datasets["elife"]
    assert elife.samples["condition"].nunique() == 11
    assert elife.samples.groupby("condition").size().eq(3).all()
    assert elife.samples["condition"].drop_duplicates().tolist() == [
        "NBF",
        "3hBF",
        "6hBF",
        "12hBF",
        "24hBF",
        "48hBF",
        "72hBF",
        "96hBF",
        "6dBF.Retained",
        "6dBF.Laid",
        "13dBF",
    ]

    crop = datasets["crop"]
    assert crop.samples["condition"].eq("NBF").all()
    assert crop.samples["tissue"].eq("crop").all()

    paper = datasets["ovary_paper"]
    assert paper.samples["reproductive_state"].nunique() == 11


def test_ir25a_resolves_in_every_dataset():
    datasets = load_datasets(EXPRESSION_DIR)
    for dataset in datasets.values():
        matches = search_genes(dataset, "ir25a", mode="exact")
        assert len(matches) == 1
        assert matches.iloc[0]["canonical_symbol"] == "Ir25a"


def test_orco_historical_aliases_resolve_in_every_dataset():
    datasets = load_datasets(EXPRESSION_DIR)
    for dataset in datasets.values():
        for query in ("ORCO", "AaegOr7", "AAEL005776"):
            matches = search_genes(dataset, query, mode="exact")
            assert len(matches) == 1
            assert matches.iloc[0]["canonical_symbol"] == "Orco"


def test_ir_family_and_expression_long():
    dataset = load_datasets(EXPRESSION_DIR)["neuro_ru"]
    members = family_members(dataset, "Ionotropic receptors (IR)")
    assert len(members) == 74
    selected = members[members["canonical_symbol"].eq("Ir25a")]
    long = expression_long(dataset, selected)
    assert len(long) == 122
    assert long["tpm"].max() > 0
    assert long["tissue"].notna().all()


def test_paper_orthology_annotations_are_available():
    datasets = load_datasets(EXPRESSION_DIR)
    for key in ("neuro_ru", "neuro_legacy"):
        orco = search_genes(datasets[key], "Orco", "exact").iloc[0]
        ir25a = search_genes(datasets[key], "Ir25a", "exact").iloc[0]
        assert orco["drosophila_ortholog"] == "Orco"
        assert ir25a["drosophila_ortholog"] == "Ir25a"


def test_nfcore_merged_gene_tpm_import(tmp_path):
    matrix = tmp_path / "salmon.merged.gene_tpm.tsv"
    matrix.write_text(
        "gene_id\tgene_name\tsample_a\tsample_b\n"
        "AAEL005776\tOrco\t12.5\t30.0\n"
        "AAEL009813\tIr25a\t2.0\t4.0\n"
    )
    dataset = load_nfcore_dataset(
        matrix,
        "nfcore_test",
        {"AAEL005776": "AaegOr7", "AAEL009813": "AaegIr25a"},
    )
    assert dataset.values.shape == (2, 2)
    assert search_genes(dataset, "Orco", "exact").iloc[0]["stable_id"] == "AAEL005776"
    assert search_genes(dataset, "Ir25a", "exact").iloc[0]["stable_id"] == "AAEL009813"
