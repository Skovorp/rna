from pathlib import Path

import pandas as pd
import pytest

from analysis.compare_elife_tpm import build_gene_map
from expression_explorer.data import load_datasets, search_genes
from expression_explorer.gene_aliases import PublishedAliases
from expression_explorer.mapping_audit import mapping_source_summary, mapping_table


EXPRESSION = Path(__file__).resolve().parents[2] / "expression"


def test_source_counts_distinguish_identifier_links_names_and_local_ids():
    summary = mapping_source_summary(PublishedAliases.load(EXPRESSION), EXPRESSION).set_index("source")
    goldman = summary.loc["goldman_2025_s1.4"]
    assert (goldman.source_rows, goldman.loc_aael_pairs, goldman.other_identifier_pairs) == (1883, 1619, 379)
    matthews = summary.loc["matthews_2016_s5"]
    assert (matthews.id_name_pairs, matthews.primary_identifiers) == (558, 554)
    morita = summary.loc["morita_2025"]
    assert (morita.id_name_pairs, morita.id_name_pairs_also_in_goldman,
            morita.id_name_pairs_only_in_this_source) == (56, 50, 5)
    assert (goldman.ppk_names, matthews.ppk_names, morita.ppk_names) == (32, 37, 31)
    jove = summary.loc["jove_2020"]
    assert (jove.identifier_pairs, jove.dataset_local_gene_loc_pairs, jove.id_name_pairs) == (0, 14412, 264)
    assert summary.loc["venkataraman_2023", "id_name_pairs"] == 0


def test_source_counts_deduplicate_rows_without_claiming_transitive_links(tmp_path):
    aliases = PublishedAliases(dict(policy="explicit_paper_rows_only", sources={"a": {}, "b": {}}, records=[
        dict(source="a", identifiers=["AAEL000001", "LOC100"], names=["ppk1"], primary="AAEL000001"),
        dict(source="a", identifiers=["loc100", "aael000001"], names=["PPK1"], primary="aael000001"),
        dict(source="b", identifiers=["AAEL000001"], names=["ppk1", "ppk2"], primary="AAEL000001"),
    ]))
    summary = mapping_source_summary(aliases, tmp_path).set_index("source")
    assert summary.loc["a", "identifier_pairs"] == 1
    assert summary.loc["a", "ppk_names"] == 1
    assert summary.id_name_pairs.tolist() == [2, 2]
    # LOC100/ppk2 is reachable in search, but neither paper asserts that pair.
    assert summary.id_name_pairs_only_in_this_source.tolist() == [1, 1]


def registry(records):
    return PublishedAliases(dict(policy="explicit_paper_rows_only", sources={}, records=[
        dict(source="test_paper", row=i + 2, primary=ids[0], identifiers=ids, names=names,
             preferred=names[0] if names else "") for i, (ids, names) in enumerate(records)
    ]))


def test_equal_symbols_and_similar_ids_do_not_create_identifier_links():
    aliases = registry([(["AAEL000001", "LOC100"], ["Shared"]),
                        (["AAEL000002"], ["Shared"]),
                        (["AAEL000001.1"], ["Different"])])
    assert aliases.identity("AAEL000001") == aliases.identity("LOC100")
    assert aliases.identity("AAEL000001") != aliases.identity("AAEL000002")
    assert aliases.identity("AAEL000001") != aliases.identity("AAEL000001.1")
    assert aliases.identity("Shared") is None
    assert aliases.identity("gene123") is None
    assert aliases.describe("AAEL000001")["mapping_status"] == "ambiguous published aliases"
    assert "AAEL000001" in aliases.describe("AAEL000001")["display_name"]


def test_comparison_uses_published_ids_and_excludes_multiple_targets():
    aliases = registry([(["LOC100", "AAEL000001"], ["PublishedName"]),
                        (["LOC200", "FeatureA", "FeatureB"], ["Split"]),
                        (["AAEL000003"], ["CommonName"])])
    published = pd.DataFrame({"IDs": ["LOC100", "LOC200", "CommonName", "Unknown", "Same"], "Symbols": [""] * 5})
    target = pd.DataFrame({"IDs": ["AAEL000001", "FeatureA", "FeatureB", "AAEL000003", "Same"], "Symbols": [""] * 5})
    mapped, diagnostics = build_gene_map(published, target, aliases)
    assert mapped[["published_id", "reanalysis_id"]].values.tolist() == [["LOC100", "AAEL000001"], ["Same", "Same"]]
    assert set(mapped.mapping_method) == {"explicit_published_identifiers", "direct_identifier"}
    assert diagnostics["ambiguous_source_aliases"] == 1
    assert 'test_paper' in mapped.iloc[0].mapping_evidence


def test_many_sources_to_one_target_are_excluded():
    aliases = registry([(["LOC1", "LOC2", "Target"], [])])
    published = pd.DataFrame({"IDs": ["LOC1", "LOC2", "Same"], "Symbols": [""] * 3})
    target = pd.DataFrame({"IDs": ["Target", "Same"], "Symbols": [""] * 2})
    mapped, diagnostics = build_gene_map(published, target, aliases)
    assert mapped.published_id.tolist() == ["Same"]
    assert diagnostics["sources_removed_for_duplicate_target"] == 2


def test_generated_mapping_policy_is_rejected():
    with pytest.raises(ValueError, match="explicit paper rows"):
        PublishedAliases(dict(policy="coordinate_matching", sources={}, records=[]))


def test_identical_paper_local_ids_do_not_match_across_studies():
    frame = pd.DataFrame({"IDs": ["gene123", "AAEL000001"], "Symbols": ["", ""]})
    mapped, diagnostics = build_gene_map(frame, frame.copy(), registry([]))
    assert mapped.published_id.tolist() == ["AAEL000001"]
    assert diagnostics["scoped_internal_ids_excluded"] == 1


def test_ppk317_resolves_by_published_aliases_with_both_descriptions():
    datasets = load_datasets(EXPRESSION)
    for dataset in datasets.values():
        for query in ("Ppk317", "ppk00873", "AAEL000873", "LOC5567199"):
            matched = search_genes(dataset, query)
            assert len(matched) == 1, (dataset.key, query)
            row = matched.iloc[0]
            assert row.display_name == "ppk317"
            assert row.vectorbase_description == "pickpocket 317"
            assert row.ncbi_description == "pickpocket protein"
            assert "goldman_2025_s1.4" in row.mapping_sources
    # Original source IDs survive, including the Jové paper's local row ID.
    assert search_genes(datasets["jove"], "Ppk317").iloc[0].stable_id == "gene16643"
    assert search_genes(datasets["ovary_paper"], "Ppk317").iloc[0].stable_id == "LOC5567199"


def test_missing_descriptions_and_internal_identifiers_are_not_guessed():
    aliases = PublishedAliases.load(EXPRESSION)
    unknown = aliases.describe("UNREPORTED_ID", "invented name")
    assert unknown["aliases"] == ()
    assert unknown["vectorbase_description"] == unknown["ncbi_description"] == ""
    assert aliases.identity("manual annotation") is None
    assert aliases.identity("gene7417") is None
    # This LOC/AAEL link existed only in the removed coordinate map.
    assert aliases.identity("LOC110674347") is None


def test_mapping_export_preserves_every_visible_row_and_protects_private_studies():
    datasets = load_datasets(EXPRESSION)
    public = mapping_table(datasets)
    all_rows = mapping_table(datasets, unlocked=True)
    assert not {"crop", "yedlin"} & set(public.dataset)
    assert len(all_rows) == sum(len(d.genes) for d in datasets.values())
    assert len(public) == sum(len(d.genes) for k, d in datasets.items() if k not in {"crop", "yedlin"})
    assert {"original_gene_id", "original_gene_name", "mapping_evidence", "searchable_aliases",
            "vectorbase_description", "ncbi_description"} <= set(public)
    ppk = public[public.display_name == "ppk317"]
    assert len(ppk) == len(datasets) - 2
    assert ppk.mapping_evidence.str.contains("goldman_2025_s1.4").all()
    assert (public.mapping_status == "original identifiers only").any()
