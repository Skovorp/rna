from pathlib import Path

from streamlit.testing.v1 import AppTest


APP = Path(__file__).resolve().parents[1] / "app.py"


def _widgets_with_options(app, expected_options):
    matches = []
    for widget_type in ("segmented_control", "radio", "pills", "selectbox"):
        for widget in getattr(app, widget_type):
            if list(widget.options) == expected_options:
                matches.append(widget)
    return matches


def _rendered_gene_names(app):
    names = set()
    for element in [*app.dataframe, *app.table]:
        frame = element.value
        for column in getattr(frame, "columns", []):
            if "gene" in str(column).casefold():
                names.update(str(value).strip().casefold() for value in frame[column])
    return names


def test_default_app_renders_without_exceptions(monkeypatch):
    monkeypatch.syspath_prepend(str(APP.parent))
    app = AppTest.from_file(str(APP), default_timeout=45).run()

    assert not app.exception, [exception.message for exception in app.exception]
    assert [title.value for title in app.title] == ["Aedes RNA Atlas"]
    assert not app.tabs

    mode_selectors = _widgets_with_options(app, ["Genes", "Families"])
    assert len(mode_selectors) == 1
    assert mode_selectors[0].value == "Genes"

    query_text = " ".join(
        widget.value
        for widget in [*app.text_input, *app.text_area]
        if isinstance(widget.value, str)
    )
    for separator in ",;/\n":
        query_text = query_text.replace(separator, " ")
    assert {"ir25a", "orco"} <= {
        token.casefold() for token in query_text.split()
    }
    assert {"ir25a", "orco"} <= _rendered_gene_names(app)

    widget_types = (
        "button",
        "checkbox",
        "download_button",
        "file_uploader",
        "multiselect",
        "number_input",
        "pills",
        "radio",
        "segmented_control",
        "select_slider",
        "selectbox",
        "slider",
        "text_area",
        "text_input",
        "toggle",
    )
    widget_count = sum(len(getattr(app, widget_type)) for widget_type in widget_types)
    table_count = len(app.dataframe) + len(app.table)
    assert widget_count + table_count <= 20


def test_public_mode_hides_persistent_import(monkeypatch):
    monkeypatch.syspath_prepend(str(APP.parent))
    monkeypatch.setenv("RNA_ATLAS_PUBLIC", "1")
    app = AppTest.from_file(str(APP), default_timeout=45).run()
    assert not app.exception
    assert all(button.label != "Import nf-core TPM" for button in app.button)
