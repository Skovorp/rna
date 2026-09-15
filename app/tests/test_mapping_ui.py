from io import BytesIO
from pathlib import Path
from zipfile import ZipFile

import pandas as pd
from streamlit.testing.v1 import AppTest


APP = Path(__file__).resolve().parents[1]


def test_methods_download_contains_gene_rows_sources_and_current_comparison_pairs(monkeypatch):
    monkeypatch.syspath_prepend(str(APP))
    monkeypatch.setattr("streamlit.page_link", lambda *args, **kwargs: None)
    downloads = []

    def capture_download(*args, **kwargs):
        downloads.append((args, kwargs))
        return False

    monkeypatch.setattr("streamlit.download_button", capture_download)
    app = AppTest.from_file(str(APP / "pages/3_Methods.py"), default_timeout=60).run()
    assert not app.exception, [item.message for item in app.exception]
    assert app.title[0].value == "Methods and gene mappings"
    body = " ".join(item.value for item in app.markdown)
    assert "Table S1.4" in body and "VectorBase description" in body and "NCBI description" in body
    assert "1,619 LOC ↔ AAEL pairs" in body and "558 ID ↔ name pairs" in body
    assert "50 also in Goldman, 1 more in Matthews, 5 only in Morita" in body
    assert "14,412" in body and "this dataset only" in body
    args, kwargs = downloads[0]
    assert kwargs["file_name"] == "aedes_gene_mapping.zip"
    with ZipFile(BytesIO(args[1])) as archive:
        assert {"gene_mapping.tsv", "published_source_rows.tsv", "mapping_source_summary.tsv",
                "processing_methods.md", "sources.json",
                "ovary_comparison_pairs.tsv", "neurotranscriptome_comparison_pairs.tsv"} <= set(archive.namelist())
        summary = pd.read_csv(archive.open("mapping_source_summary.tsv"), sep="\t").set_index("source")
        assert summary.loc["goldman_2025_s1.4", "loc_aael_pairs"] == 1619
        assert summary.loc["jove_2020", "dataset_local_gene_loc_pairs"] == 14412
        assert "max_memory" in archive.read("processing_methods.md").decode()
        genes = pd.read_csv(archive.open("gene_mapping.tsv"), sep="\t", keep_default_na=False)
        assert not {"yedlin", "crop"} & set(genes.dataset)
        assert genes.original_gene_id.ne("").all()
        assert genes.matrix_file.ne("").all()
        assert genes.loc[genes.display_name.eq("ppk317"), "ncbi_description"].eq("pickpocket protein").all()
        for study in ("ovary", "neurotranscriptome"):
            pairs = pd.read_csv(archive.open(f"{study}_comparison_pairs.tsv"), sep="\t")
            assert set(pairs.mapping_method) <= {"direct_identifier", "explicit_published_identifiers"}
            assert pairs.mapping_evidence.notna().all()


def test_gene_details_show_both_original_descriptions():
    app = AppTest.from_file(str(APP / "app.py"), default_timeout=60)
    app.query_params["page"] = "Genes"
    app.run()
    query = next(widget for widget in app.text_input if widget.label == "Genes or identifiers")
    query.set_value("Ppk317").run()
    assert not app.exception, [item.message for item in app.exception]
    descriptions = next(expander for expander in app.expander if expander.label == "ppk317")
    text = " ".join(item.value for item in descriptions.markdown)
    assert "VectorBase" in text and "pickpocket 317" in text
    assert "NCBI" in text and "pickpocket protein" in text
    captions = " ".join(item.value for item in descriptions.caption)
    assert "AAEL000873" in captions and "LOC5567199" in captions
