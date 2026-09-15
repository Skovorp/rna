from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

from expression_explorer.catalog import tpm_download_path, visible_catalog
from expression_explorer.data import load_datasets, search_genes


APP = Path(__file__).resolve().parents[1] / "app.py"
EXPRESSION = APP.parents[1] / "expression"


def test_private_catalog_and_downloads_require_unlock():
    assert not {"crop", "yedlin"} & {entry.key for entry in visible_catalog()}
    assert {"crop", "yedlin"} <= {entry.key for entry in visible_catalog(True)}
    for key in ("crop", "yedlin"):
        with pytest.raises(PermissionError):
            tpm_download_path(key, EXPRESSION)
        assert tpm_download_path(key, EXPRESSION, unlocked=True).is_file()
    for key in ("elife", "atlas", "midgut"):
        assert tpm_download_path(key, EXPRESSION).is_file()
    with pytest.raises(KeyError):
        tpm_download_path("../METHODS.md", EXPRESSION)


def test_home_downloads_and_private_prompt(monkeypatch):
    monkeypatch.syspath_prepend(str(APP.parent))
    app = AppTest.from_file(str(APP), default_timeout=45).run()
    assert not app.exception
    downloads = app.get("download_button")
    assert {item.key for item in downloads} == {"download_tpm_elife", "download_tpm_atlas", "download_tpm_midgut"}
    html = " ".join(element.value for element in app.markdown)
    assert html.index("## Datasets") < html.index("Fruitless exons (published)") < html.index("Methodology →")
    assert all(accession in html for accession in ("PRJNA796320", "PRJNA236239", "PRJNA1020561", "PRJNA605870", "PRJNA612100"))
    app.button(key="catalog_unlock_private").click().run()
    assert not app.exception
    assert app.text_input[0].label == "Password"
    app.text_input[0].set_value("incorrect-password-for-regression-test").run()
    assert any(item.value == "Wrong password." for item in app.error)
    assert "private_datasets_unlocked" not in app.session_state or not app.session_state["private_datasets_unlocked"]
    app.session_state["private_datasets_unlocked"] = True
    app.session_state["site_navigation"] = "Home"
    app.run()
    assert not app.exception
    assert {item.key for item in app.get("download_button")} == {
        "download_tpm_elife", "download_tpm_atlas", "download_tpm_midgut",
        "download_tpm_crop", "download_tpm_yedlin",
    }
    assert not any(button.key == "catalog_unlock_private" for button in app.button)


def test_published_tpm_dimensions_units_and_aliases():
    datasets = load_datasets(EXPRESSION)
    for key, year, shape in (("morita", 2025, (18943, 10)), ("jove", 2020, (14676, 12))):
        dataset = datasets[key]
        assert dataset.values.shape == shape
        source = pd.read_csv(EXPRESSION / f"{key}_{year}_gene_tpm.tsv.gz", sep="\t")
        np.testing.assert_array_equal(dataset.values, source[dataset.sample_columns])
        for query in ("Ir25a", "Orco", "AAEL005776", "Ir7a"):
            assert len(search_genes(dataset, query)) == 1
        assert dataset.samples.groupby("condition_label").size().nunique() == 1
    assert datasets["morita"].samples.genotype.nunique() == 2
    assert datasets["jove"].samples.tissue.unique().tolist() == ["stylet", "labium"]


@pytest.mark.parametrize("key", ["morita", "jove"])
def test_new_study_deep_link_opens_gene_explorer(monkeypatch, key):
    monkeypatch.syspath_prepend(str(APP.parent))
    app = AppTest.from_file(str(APP), default_timeout=45)
    app.query_params.update({"page": "Genes", "study": key})
    app.run()
    assert not app.exception, [item.message for item in app.exception]
    assert app.multiselect(key="gene_studies").value == [key]
    assert app.get("plotly_chart")


@pytest.mark.parametrize("key", ["morita", "jove"])
@pytest.mark.parametrize("page", ["Families", "Clusters"])
def test_published_studies_work_in_other_tpm_explorers(monkeypatch, key, page):
    monkeypatch.syspath_prepend(str(APP.parent))
    app = AppTest.from_file(str(APP), default_timeout=45)
    app.query_params["page"] = page
    if page == "Families":
        app.session_state["family_studies"] = [key]
    else:
        app.session_state["cluster_study"] = key
    app.run()
    assert not app.exception, [item.message for item in app.exception]
    assert app.get("plotly_chart")


def test_basrur_preserves_missing_replicates_and_exon_units(monkeypatch):
    monkeypatch.syspath_prepend(str(APP.parent))
    data = pd.read_csv(EXPRESSION / "basrur_2020_fruitless_exon_counts.tsv.gz", sep="\t")
    assert data.groupby("panel").size().to_dict() == {"Figure 1G": 273, "Figure 1H": 144}
    brain = data[data.panel.eq("Figure 1G")]
    assert brain.groupby(["exon", "sex"]).size().unstack()[["female", "male"]].eq([4, 3]).all().all()
    first = brain[(brain.exon_number == 1) & brain.sex.eq("female") & brain.replicate.eq(1)]
    assert first.normalized_count.iloc[0] == pytest.approx(153.846153846154)
    assert not (data.tissue.str.contains("ovar", case=False) & data.sex.eq("male")).any()
    app = AppTest.from_file(str(APP), default_timeout=45)
    app.query_params.update({"page": "Genes", "view": "fruitless"})
    app.run()
    assert not app.exception, [item.message for item in app.exception]
    assert len(app.get("plotly_chart")) == 1
    assert not any(widget.label == "TPM scale" for widget in app.selectbox)
    app.radio[0].set_value("Brain (all exons)").run()
    assert not app.exception
    assert len(app.dataframe[0].value) == 273
