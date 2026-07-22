import json
from pathlib import Path
from statistics import median

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


def _plotly_spec(app, index=0):
    return json.loads(app.get("plotly_chart")[index].proto.spec)


def test_default_app_renders_without_exceptions(monkeypatch):
    monkeypatch.syspath_prepend(str(APP.parent))
    app = AppTest.from_file(str(APP), default_timeout=45).run()

    assert not app.exception, [exception.message for exception in app.exception]
    assert [title.value for title in app.title] == ["Aedes RNA Atlas"]
    assert not app.tabs
    assert 'label="Mosquito basics"' in APP.read_text()

    mode_selectors = _widgets_with_options(
        app, ["Genes", "Families", "Compare conditions"]
    )
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
    assert any(
        "Samples ≥1 TPM" in element.value.columns
        for element in app.dataframe
    )

    studies = next(widget for widget in app.multiselect if widget.label == "Studies")
    assert len(studies.options) == 2
    assert all("legacy" not in option.casefold() for option in studies.options)
    assert studies.value == ["elife", "neuro_ru"]

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


def test_gene_plots_are_horizontal_and_sortable(monkeypatch):
    monkeypatch.syspath_prepend(str(APP.parent))
    app = AppTest.from_file(str(APP), default_timeout=45).run()
    assert not app.exception

    plot = _plotly_spec(app)
    assert plot["data"][0]["orientation"] == "h"
    assert plot["layout"]["xaxis"]["title"]["text"] == "log₂(TPM + 1)"
    assert plot["layout"]["xaxis"]["autorange"] is False
    assert plot["layout"]["xaxis"]["range"][0] == 0
    assert plot["layout"]["xaxis"]["range"][1] > 0
    assert plot["layout"]["yaxis"]["title"]["text"] == "Reproductive state"
    assert plot["layout"]["yaxis"]["autorange"] == "reversed"
    assert plot["layout"]["legend"]["itemclick"] is False
    assert plot["layout"]["legend"]["itemdoubleclick"] is False
    median_trace = next(trace for trace in plot["data"] if trace.get("name") == "Group median")
    assert median_trace["showlegend"] is False

    condition_order = plot["layout"]["yaxis"]["categoryarray"]
    guides = plot["layout"]["shapes"]
    assert [guide["y0"] for guide in guides] == condition_order[::2]
    assert all(guide["line"]["dash"] == "dot" for guide in guides)
    assert all(guide["layer"] == "below" for guide in guides)

    median_toggle = next(
        toggle for toggle in app.toggle if toggle.label == "Show group medians"
    )
    median_toggle.set_value(False).run()
    without_medians = _plotly_spec(app)
    assert all(trace.get("name") != "Group median" for trace in without_medians["data"])
    assert without_medians["layout"]["height"] == plot["layout"]["height"]
    assert without_medians["layout"]["xaxis"]["range"] == plot["layout"]["xaxis"]["range"]
    assert (
        without_medians["layout"]["yaxis"]["categoryarray"]
        == plot["layout"]["yaxis"]["categoryarray"]
    )

    guide_toggle = next(
        toggle for toggle in app.toggle if toggle.label == "Show row guides"
    )
    guide_toggle.set_value(False).run()
    without_guides = _plotly_spec(app)
    assert without_guides["layout"].get("shapes", []) == []
    assert without_guides["layout"]["height"] == plot["layout"]["height"]
    assert without_guides["layout"]["xaxis"]["range"] == plot["layout"]["xaxis"]["range"]
    assert (
        without_guides["layout"]["yaxis"]["categoryarray"]
        == plot["layout"]["yaxis"]["categoryarray"]
    )

    sort_toggle = next(
        toggle for toggle in app.toggle if toggle.label == "Sort conditions by expression"
    )
    sort_toggle.set_value(True).run()
    sorted_plot = _plotly_spec(app)
    first_trace = sorted_plot["data"][0]
    values_by_condition = {}
    for condition, customdata in zip(first_trace["y"], first_trace["customdata"]):
        values_by_condition.setdefault(condition, []).append(customdata[1])
    expected_order = sorted(
        values_by_condition,
        key=lambda condition: median(values_by_condition[condition]),
        reverse=True,
    )
    assert sorted_plot["layout"]["yaxis"]["categoryarray"] == expected_order

    studies = next(widget for widget in app.multiselect if widget.label == "Studies")
    studies.set_value(["elife", "neuro_ru"]).run()
    assert not app.exception
    study_plots = [
        _plotly_spec(app, index)
        for index in range(len(app.get("plotly_chart")))
        if _plotly_spec(app, index)["data"][0]["type"] == "box"
    ]
    assert len(study_plots) == 4
    assert all(plot["data"][0]["orientation"] == "h" for plot in study_plots)
    comparison_heatmaps = [
        _plotly_spec(app, index)
        for index in range(len(app.get("plotly_chart")))
        if _plotly_spec(app, index)["data"][0]["type"] == "heatmap"
    ]
    assert len(comparison_heatmaps) == 2


def test_family_mode_explains_filter_and_zscore(monkeypatch):
    monkeypatch.syspath_prepend(str(APP.parent))
    app = AppTest.from_file(str(APP), default_timeout=45).run()
    mode = _widgets_with_options(
        app, ["Genes", "Families", "Compare conditions"]
    )[0]
    mode.set_value("Families").run()

    assert not app.exception
    captions = " ".join(element.value for element in app.caption)
    assert "does not combine the family into one score" in captions
    zscore_toggle = next(
        toggle
        for toggle in app.toggle
        if toggle.label == "Show relative pattern within each gene"
    )
    assert "z-scores" in zscore_toggle.help


def test_condition_comparison_uses_fdr(monkeypatch):
    monkeypatch.syspath_prepend(str(APP.parent))
    app = AppTest.from_file(str(APP), default_timeout=45).run()
    mode = _widgets_with_options(
        app, ["Genes", "Families", "Compare conditions"]
    )[0]
    mode.set_value("Compare conditions").run()

    assert not app.exception, [exception.message for exception in app.exception]
    result_tables = [frame.value for frame in app.dataframe if "FDR" in frame.value.columns]
    assert len(result_tables) == 1
    assert {
        "Gene",
        "Mean TPM (A)",
        "Mean TPM (B)",
        "Average TPM",
        "TPM ratio (A / B)",
        "Raw p-value",
        "FDR",
        "FDR < 0.05",
    } <= set(result_tables[0].columns)
    assert len(result_tables[0]) > 10_000

    plot = _plotly_spec(app)
    assert len(plot["data"]) == 1
    assert plot["data"][0]["type"] == "scattergl"
    assert plot["data"][0]["showlegend"] is False
    assert plot["data"][0]["marker"] == {
        "color": "#f5b85b",
        "opacity": 1.0,
        "size": 3,
    }
    assert plot["layout"]["xaxis"]["title"]["text"] == "Average TPM (logarithmic scale)"
    assert plot["layout"]["xaxis"]["type"] == "log"
    assert plot["layout"]["xaxis"]["range"][0] < plot["layout"]["xaxis"]["range"][1]
    assert {"1", "10", "100", "1,000"} <= set(plot["layout"]["xaxis"]["ticktext"])
    assert (
        plot["layout"]["yaxis"]["title"]["text"]
        == "TPM ratio A / B (logarithmic scale)"
    )
    assert plot["layout"]["yaxis"]["type"] == "log"
    assert plot["layout"]["yaxis"]["range"][0] == -plot["layout"]["yaxis"]["range"][1]
    assert {"0.1×", "1×", "10×"} <= set(plot["layout"]["yaxis"]["ticktext"])
    assert plot["layout"]["shapes"][0]["y0"] == 1
    assert plot["layout"]["shapes"][0]["y1"] == 1
    captions = " ".join(element.value for element in app.caption)
    assert "A/B ratio is undefined" in captions
    assert "fully opaque and uses the same color regardless of FDR" in captions
    fdr_threshold = next(
        widget for widget in app.number_input if widget.label == "FDR threshold"
    )
    assert fdr_threshold.value == 0.05
    fdr_threshold.set_value(0.1).run()
    assert any("FDR < 0.1" in frame.value.columns for frame in app.dataframe)
    assert any(button.label == "Download all comparison results" for button in app.download_button)
