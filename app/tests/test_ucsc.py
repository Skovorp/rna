from pathlib import Path

from expression_explorer.ucsc import (
    cell_browser_expression_url,
    cell_browser_metadata_url,
    cell_browser_url,
    find_gene_matches,
    load_manifest,
)


MANIFEST = (
    Path(__file__).resolve().parents[2]
    / "expression"
    / "ucsc_mosquito_cell_atlas_genes.json.gz"
)


def test_cell_browser_url_matches_ucsc_share_format():
    assert cell_browser_url("mosquito/all", "Orco") == (
        "/ucsc/?ds=mosquito+all&gene=Orco"
    )


def test_cell_browser_expression_url_supports_multiple_genes_and_grouping():
    assert cell_browser_expression_url(
        "mosquito/all",
        ["Ir25a", "Orco", "AAEL021429"],
        "annotation",
        context_gene="Ir25a",
    ) == (
        "/ucsc/?ds=mosquito+all&gene=Ir25a&"
        "exprGene=Ir25a+Orco+AAEL021429&exprMeta=annotation"
    )


def test_cell_browser_metadata_url_matches_ucsc_share_format():
    assert cell_browser_metadata_url("mosquito/t012", "sample") == (
        "/ucsc/?ds=mosquito+t012&meta=sample"
    )


def test_manifest_resolves_symbols_when_legacy_ids_are_not_indexed():
    datasets = load_manifest(MANIFEST)
    assert len(datasets) == 24
    assert datasets[0].name == "mosquito/all"
    assert datasets[0].sample_count == 330_364
    assert datasets[0].default_metadata_field == "annotation"
    assert ("annotation", "annotation") in datasets[0].categorical_fields
    assert datasets[0].quick_genes[0] == ("AAEL021429|spir", "Unspecified")

    orco = find_gene_matches(datasets, ["AAEL005776", "AaegOr7", "Orco"])
    ir25a = find_gene_matches(datasets, ["AAEL009813", "AaegIr25a", "Ir25a"])

    assert len(orco) == 24
    assert len(ir25a) == 24
    assert orco[0].gene_query == "Orco"
    assert ir25a[0].gene_query == "Ir25a"
