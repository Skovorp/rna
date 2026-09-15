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
    assert not app.get("download_button")
    html = " ".join(element.value for element in app.markdown)
    assert html.count('>Download TPM tables</a>') == 3
    assert 'download="crop_star_salmon_gene_tpm.tsv.gz"' not in html
    assert 'download="yedlin_star_salmon_gene_tpm.tsv.gz"' not in html
    assert html.index("## Datasets") < html.index("Mouthparts (published)") < html.index("Methodology →")
    assert "Basrur et al. (2020)" in html
    assert "Fruitless exons" not in html
    assert all(accession in html for accession in ("PRJNA796320", "PRJNA236239", "PRJNA1020561", "PRJNA605870"))
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
    html = " ".join(element.value for element in app.markdown)
    assert not app.get("download_button")
    assert html.count('>Download TPM tables</a>') == 5
    assert 'download="crop_star_salmon_gene_tpm.tsv.gz"' in html
    assert 'download="yedlin_star_salmon_gene_tpm.tsv.gz"' in html
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


def test_old_basrur_link_opens_reprocessed_gene_tpm(monkeypatch):
    monkeypatch.syspath_prepend(str(APP.parent))
    app = AppTest.from_file(str(APP), default_timeout=45)
    app.query_params.update({"page": "Genes", "view": "fruitless"})
    app.run()
    assert not app.exception, [item.message for item in app.exception]
    assert app.multiselect(key="gene_studies").value == ["atlas"]
    assert any(widget.label == "TPM scale" for widget in app.selectbox)
    assert "Fruitless exon" not in " ".join(element.value for element in app.markdown)
    assert app.get("plotly_chart")
