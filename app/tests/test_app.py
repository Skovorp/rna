import json
from pathlib import Path

import numpy as np
from streamlit.testing.v1 import AppTest


APP = Path(__file__).resolve().parents[1] / "app.py"
NAVIGATION_ITEMS = ["Home", "Genes", "Families", "Differential expression", "Clusters"]
MATCHED_GENE_COLUMNS = {
    "Include",
    "Gene",
    "Family",
    "Alternative names",
    "Study coverage",
}


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


def _plotly_specs_by_type(app, chart_type):
    return [
        _plotly_spec(app, index)
        for index in range(len(app.get("plotly_chart")))
        if _plotly_spec(app, index)["data"][0]["type"] == chart_type
    ]


def _scatter_genes(plot):
    return {
        str(trace["name"])
        for trace in plot["data"]
        if trace["type"] == "box"
    }


def _matched_gene_editor(app):
    matches = [
        element
        for element in app.dataframe
        if MATCHED_GENE_COLUMNS <= set(element.value.columns)
        and any(
            str(column).startswith("Mean TPM: ")
            for column in element.value.columns
        )
    ]
    assert len(matches) == 1
    return matches[0]


def _mean_tpm_columns(frame):
    return [
        column for column in frame.columns if str(column).startswith("Mean TPM: ")
    ]


def _family_gene_tables(app):
    return [
        element.value
        for element in app.dataframe
        if {"Gene", "Stable ID", "Mean TPM", "Top context"}
        <= set(element.value.columns)
    ]


def _set_matched_gene_enabled(app, gene, enabled):
    editor = _matched_gene_editor(app)
    row_index = editor.value.index[editor.value["Gene"] == gene].item()
    app.session_state[editor.key] = {
        "edited_rows": {int(row_index): {"Include": enabled}},
        "deleted_rows": [],
        "added_rows": [],
    }
    return app.run()


def _select_page(app, page):
    navigation = _widgets_with_options(app, NAVIGATION_ITEMS)
    assert len(navigation) == 1
    navigation[0].set_value(page).run()
    return app


def test_default_app_renders_without_exceptions(monkeypatch):
    monkeypatch.syspath_prepend(str(APP.parent))
    app = AppTest.from_file(str(APP), default_timeout=45).run()

    assert not app.exception, [exception.message for exception in app.exception]
    assert not app.title
    assert not app.tabs

    mode_selectors = _widgets_with_options(app, NAVIGATION_ITEMS)
    assert len(mode_selectors) == 1
    assert mode_selectors[0].value == "Home"
    assert app.query_params["page"] == ["Home"]
    home_html = " ".join(element.value for element in app.markdown)
    assert "# Aedes RNA Atlas" in home_html
    assert "Use the menu above to:" in home_html
    assert "## Datasets" in home_html
    assert "## Comparisons" in home_html
    assert "**Ovary (paper)**" in home_html
    assert "**Ovary (reprocessed)**" in home_html
    assert "**Atlas (paper)**" in home_html
    assert "**Atlas (reprocessed)**" in home_html
    assert "**Midgut (reprocessed)**" in home_html
    assert "**Crop (reprocessed)**" in home_html
    assert "identical* pipeline" in home_html
    assert "never recomputed from TPM" in home_html
    # Comparisons and Methods are reachable as inline links, not just from the
    # page-link buttons at the bottom.
    assert home_html.count("(/Ovary_paper_vs_reprocessed)") >= 3
    assert home_html.count("(/Atlas_paper_vs_reprocessed)") >= 3
    assert "(/Methods)" in home_html
    # The bottom page-link duplicate of the ovary comparison was removed.
    assert "Ovary: paper vs reprocessed comparison" not in APP.read_text()
    # Pipeline parameters live only in METHODS.md, rendered by the Methods page.
    assert "nextflow run" not in home_html
    assert "skip_alignment" not in home_html
    assert 'class="home-hero"' not in home_html
    assert 'class="home-card-grid"' not in home_html
    assert 'class="home-metrics"' not in home_html

    source = APP.read_text()
    assert "nextflow run" not in source
    assert "pseudo_aligner" not in source
    assert '[data-testid="stSidebar"]' in source
    assert '[data-testid="stSidebarCollapsedControl"]' in source
    assert '[data-testid="stToolbar"]' in source
    assert "with st.sidebar" not in source
    assert "not a third study" not in source

    _select_page(app, "Genes")
    assert not app.exception, [exception.message for exception in app.exception]
    assert app.query_params["page"] == ["Genes"]

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
    token_html = next(
        element.value
        for element in app.markdown
        if '<div class="gene-token-overlay"' in element.value
    )
    assert "Ir25a" in token_html
    assert "Orco" in token_html
    assert token_html.count("gene-token-found") == 2
    assert {"ir25a", "orco"} <= _rendered_gene_names(app)
    assert any(
        "Samples ≥1 TPM" in element.value.columns
        for element in app.dataframe
    )
    summary_tables = [
        element.value
        for element in app.dataframe
        if "Samples ≥1 TPM" in element.value.columns
    ]
    assert len(summary_tables) == 1
    assert {"ir25a", "orco"} <= {
        str(value).casefold() for value in summary_tables[0]["Gene"]
    }
    raw_tables = [
        element.value
        for element in app.dataframe
        if {"Gene", "Sample", "TPM"} <= set(element.value.columns)
    ]
    assert len(raw_tables) == 1
    assert ".st-key-gene_setup_panel" in source
    matched_genes = _matched_gene_editor(app).value
    assert matched_genes["Include"].all()
    assert set(matched_genes["Gene"]) == {"Ir25a", "Orco"}
    mean_tpm_columns = _mean_tpm_columns(matched_genes)
    assert len(mean_tpm_columns) == 2
    assert matched_genes[mean_tpm_columns[0]].dropna().is_monotonic_decreasing
    for mean_tpm_column in mean_tpm_columns:
        study_label = mean_tpm_column.removeprefix("Mean TPM: ")
        study_values = raw_tables[0][raw_tables[0]["Study"] == study_label]
        raw_mean_tpm = study_values.groupby("Gene")["TPM"].mean()
        for _, matched_gene in matched_genes.iterrows():
            gene = matched_gene["Gene"]
            if gene in raw_mean_tpm:
                assert abs(matched_gene[mean_tpm_column] - raw_mean_tpm[gene]) < 1e-12
            else:
                assert matched_gene[mean_tpm_column] != matched_gene[mean_tpm_column]
    assert any(button.label == "Turn all on" for button in app.button)
    assert any(button.label == "Turn all off" for button in app.button)
    info_html = " ".join(
        element.value
        for element in app.markdown
        if "section-title-with-info" in element.value
    )
    assert "Expression across selected genes" in info_html
    assert "sample filters" in info_html
    assert "optional condition-label coloring" in info_html
    assert "Point colors identify genes" in info_html
    assert "its heatmap" in info_html
    assert "Expression heatmaps" not in info_html
    assert 'class="section-info-tooltip"' in info_html
    assert 'tabindex="0"' in info_html
    captions = " ".join(element.value for element in app.caption)
    assert "One replicate plot per study" not in captions
    assert "One heatmap per study" not in captions
    assert "Group median TPM" not in captions

    studies = next(widget for widget in app.multiselect if widget.label == "Studies")
    assert len(studies.options) == 5
    assert "Midgut (reprocessed), blood-meal time course" in studies.options
    assert all("legacy" not in option.casefold() for option in studies.options)
    assert studies.value == ["elife", "neuro_ru"]

    filter_widgets = [
        widget for widget in app.multiselect if str(widget.key).startswith("gene_filter_")
    ]
    assert len(filter_widgets) == 4
    assert all(widget.value == [] for widget in filter_widgets)
    color_widgets = [
        widget
        for widget in app.selectbox
        if str(widget.key).startswith("gene_label_color_")
    ]
    assert len(color_widgets) == 2
    assert all(widget.value is None for widget in color_widgets)
    atlas_gene = next(widget for widget in app.selectbox if widget.label == "Atlas gene")
    atlas_view = next(widget for widget in app.selectbox if widget.label == "Atlas view")
    split_umap = next(
        widget for widget in app.selectbox if widget.label == "Additional UMAP"
    )
    expression_view = next(
        widget for widget in app.selectbox if widget.label == "Expression atlas view"
    )
    expression_grouping = next(
        widget for widget in app.selectbox if widget.label == "Group cells by"
    )
    include_dataset_genes = next(
        widget for widget in app.toggle if widget.label == "Include dataset genes"
    )
    assert atlas_gene.options == ["Ir25a", "Orco"]
    assert atlas_view.value == "mosquito/all"
    assert len(atlas_view.options) == 24
    assert split_umap.options == ["None", "sex", "sample"]
    assert split_umap.value == "None"
    assert expression_view.value == "mosquito/all"
    assert len(expression_view.options) == 24
    assert expression_grouping.value == "annotation"
    assert include_dataset_genes.value is False
    assert any(
        "/ucsc/?ds=mosquito+all&amp;gene=Ir25a" in element.value
        or "/ucsc/?ds=mosquito+all&gene=Ir25a" in element.value
        for element in app.markdown
    )
    assert any(
        "exprGene=Ir25a+Orco&amp;exprMeta=annotation" in element.value
        or "exprGene=Ir25a+Orco&exprMeta=annotation" in element.value
        for element in app.markdown
    )

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
    # UMAP and expression controls stay focused despite embedding two atlas views.
    # +1 for the private-datasets unlock password input.
    assert widget_count + table_count <= 27

    logo = next(button for button in app.button if button.label == "🧬 Aedes RNA Atlas")
    logo.click().run()
    assert _widgets_with_options(app, NAVIGATION_ITEMS)[0].value == "Home"
    assert any("# Aedes RNA Atlas" in element.value for element in app.markdown)


def test_url_page_state_restores_genes_after_a_fresh_session(monkeypatch):
    monkeypatch.syspath_prepend(str(APP.parent))
    app = AppTest.from_file(str(APP), default_timeout=45)
    app.query_params["page"] = "Genes"
    app.run()

    assert not app.exception
    navigation = _widgets_with_options(app, NAVIGATION_ITEMS)
    assert len(navigation) == 1
    assert navigation[0].value == "Genes"


def test_gene_atlas_expression_controls_update_the_embed(monkeypatch):
    monkeypatch.syspath_prepend(str(APP.parent))
    app = AppTest.from_file(str(APP), default_timeout=45).run()
    _select_page(app, "Genes")

    grouping = next(
        widget for widget in app.selectbox if widget.label == "Group cells by"
    )
    assert "tissues" in grouping.options
    grouping.set_value("tissues").run()

    include_dataset_genes = next(
        widget for widget in app.toggle if widget.label == "Include dataset genes"
    )
    include_dataset_genes.set_value(True).run()

    expression_links = [
        element.value
        for element in app.markdown
        if "Open this expression plot in a new tab" in element.value
    ]
    assert len(expression_links) == 1
    assert "exprGene=Ir25a+Orco+AAEL021429" in expression_links[0]
    assert "exprMeta=tissues" in expression_links[0]
    captions = " ".join(element.value for element in app.caption)
    assert "Showing 2 selected genes plus 15 dataset genes" in captions


def test_gene_atlas_umap_can_add_sample_metadata_plot_below(monkeypatch):
    monkeypatch.syspath_prepend(str(APP.parent))
    app = AppTest.from_file(str(APP), default_timeout=45).run()
    _select_page(app, "Genes")

    atlas_view = next(widget for widget in app.selectbox if widget.label == "Atlas view")
    atlas_view.set_value("mosquito/t012").run()

    split_umap = next(
        widget for widget in app.selectbox if widget.label == "Additional UMAP"
    )
    assert split_umap.options == ["None", "sex", "sample"]
    split_umap.set_value("sample").run()

    split_links = [
        element.value
        for element in app.markdown
        if "Open sample UMAP in a new tab" in element.value
    ]
    assert split_links == [
        "[Open sample UMAP in a new tab](/ucsc/?ds=mosquito+t012&meta=sample)"
    ]
    captions = " ".join(element.value for element in app.caption)
    assert "Expression: Ir25a" in captions
    assert "Colored by sample" in captions
    assert "primary_column, split_column = st.columns(2)" not in APP.read_text()


def test_home_has_no_import_controls(monkeypatch):
    monkeypatch.syspath_prepend(str(APP.parent))
    app = AppTest.from_file(str(APP), default_timeout=45).run()
    assert not app.exception
    assert not app.file_uploader
    assert all("import" not in button.label.casefold() for button in app.button)
    source = APP.read_text()
    assert "import_nfcore_dialog" not in source
    assert "Import local nf-core TPM" not in source


def test_gene_plots_are_horizontal_and_scaled(monkeypatch):
    monkeypatch.syspath_prepend(str(APP.parent))
    app = AppTest.from_file(str(APP), default_timeout=45).run()
    assert not app.exception
    _select_page(app, "Genes")

    plot = _plotly_spec(app)
    assert plot["data"][0]["orientation"] == "h"
    axis_scale = next(
        selectbox for selectbox in app.selectbox if selectbox.label == "TPM scale"
    )
    assert list(axis_scale.options) == ["Linear", "Log base 2", "Log base 10"]
    assert axis_scale.value == "Log base 2"
    assert plot["layout"]["xaxis"]["title"]["text"] == "TPM (log₂ scale)"
    assert plot["layout"]["xaxis"]["autorange"] is False
    assert plot["layout"]["xaxis"]["range"][0] == 0
    assert plot["layout"]["xaxis"]["range"][1] > 0
    assert {"0", "1", "2", "4"} <= set(plot["layout"]["xaxis"]["ticktext"])
    assert "TPM=" in plot["data"][0]["hovertemplate"]
    assert plot["layout"]["yaxis"]["title"]["text"] == "Reproductive state"
    assert plot["layout"]["yaxis"]["autorange"] == "reversed"
    assert plot["layout"]["legend"]["itemclick"] is False
    assert plot["layout"]["legend"]["itemdoubleclick"] is False
    median_trace = next(trace for trace in plot["data"] if trace.get("name") == "Group median")
    assert median_trace["showlegend"] is False

    condition_order = plot["layout"]["yaxis"]["categoryarray"]
    guides = plot["layout"]["shapes"]
    assert [guide["y0"] for guide in guides] == condition_order
    assert all(guide["line"]["dash"] == "dot" for guide in guides)
    assert all(guide["layer"] == "below" for guide in guides)
    gene_trace_colors = {
        trace["marker"]["color"]
        for trace in plot["data"]
        if trace["type"] == "box"
    }
    assert _scatter_genes(plot) == {"Ir25a", "Orco"}
    assert len(gene_trace_colors) == 2

    axis_scale.set_value("Log base 10").run()
    base10_plot = _plotly_specs_by_type(app, "box")[1]
    assert base10_plot["layout"]["xaxis"]["title"]["text"] == "TPM (log₁₀ scale)"
    assert {"0", "1", "10", "100"} <= set(
        base10_plot["layout"]["xaxis"]["ticktext"]
    )
    base10_heatmaps = _plotly_specs_by_type(app, "heatmap")
    assert base10_heatmaps
    assert all(
        heatmap["data"][0]["colorbar"]["title"]["text"]
        == "Median TPM<br>(log base 10 scale)"
        for heatmap in base10_heatmaps
    )
    assert any(
        {"1", "10", "100"}
        <= set(heatmap["data"][0]["colorbar"]["ticktext"])
        for heatmap in base10_heatmaps
    )
    assert all(
        "Median TPM: %{customdata:.3f}" in heatmap["data"][0]["hovertemplate"]
        for heatmap in base10_heatmaps
    )

    axis_scale.set_value("Linear").run()
    linear_plot = _plotly_spec(app)
    assert linear_plot["layout"]["xaxis"]["title"]["text"] == "TPM (linear scale)"
    assert "ticktext" not in linear_plot["layout"]["xaxis"]
    linear_heatmaps = _plotly_specs_by_type(app, "heatmap")
    assert all(
        heatmap["data"][0]["colorbar"]["title"]["text"]
        == "Median TPM<br>(linear scale)"
        for heatmap in linear_heatmaps
    )

    axis_scale.set_value("Log base 2").run()
    plot = _plotly_spec(app)

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

    assert all(
        toggle.label != "Sort conditions by expression" for toggle in app.toggle
    )

    studies = next(widget for widget in app.multiselect if widget.label == "Studies")
    studies.set_value(["elife", "neuro_ru"]).run()
    assert not app.exception
    study_plots = _plotly_specs_by_type(app, "box")
    assert len(study_plots) == 2
    assert all(plot["data"][0]["orientation"] == "h" for plot in study_plots)
    comparison_heatmaps = _plotly_specs_by_type(app, "heatmap")
    assert len(comparison_heatmaps) == 2
    assert [
        _plotly_spec(app, index)["data"][0]["type"]
        for index in range(len(app.get("plotly_chart")))
    ] == ["box", "heatmap", "box", "heatmap"]

    query = next(widget for widget in app.text_input if widget.label == "Genes or identifiers")
    query.set_value("Ir25a Orco Ir7a").run()
    assert not app.exception
    three_gene_plots = _plotly_specs_by_type(app, "box")
    assert len(three_gene_plots) == 2
    assert _scatter_genes(three_gene_plots[0]) == {"Ir25a", "Orco", "Ir7a"}
    assert _scatter_genes(three_gene_plots[1]) == {"Ir25a", "Orco"}


def test_gene_graph_filters_and_group_colors_do_not_change_details(monkeypatch):
    monkeypatch.syspath_prepend(str(APP.parent))
    app = AppTest.from_file(str(APP), default_timeout=45).run()
    _select_page(app, "Genes")

    raw_before = next(
        frame.value
        for frame in app.dataframe
        if {"Gene", "Sample", "TPM"} <= set(frame.value.columns)
    )
    sample_counts_before = raw_before.groupby("Study")["Sample"].nunique().to_dict()

    ovary_filter = next(
        widget
        for widget in app.multiselect
        if widget.key == "gene_filter_elife_reproductive_state"
    )
    ovary_filter.set_value(["Non-blood-fed"]).run()
    assert not app.exception

    chart_types = [
        _plotly_spec(app, index)["data"][0]["type"]
        for index in range(len(app.get("plotly_chart")))
    ]
    assert chart_types == ["box", "heatmap", "box", "heatmap"]
    ovary_plot = _plotly_specs_by_type(app, "box")[0]
    ovary_heatmap = _plotly_specs_by_type(app, "heatmap")[0]
    assert ovary_plot["layout"]["yaxis"]["categoryarray"] == ["Non-blood-fed"]
    assert ovary_heatmap["data"][0]["x"] == ["Non-blood-fed"]

    raw_after = next(
        frame.value
        for frame in app.dataframe
        if {"Gene", "Sample", "TPM"} <= set(frame.value.columns)
    )
    assert raw_after.groupby("Study")["Sample"].nunique().to_dict() == sample_counts_before
    summary = next(
        frame.value
        for frame in app.dataframe
        if "Samples ≥1 TPM" in frame.value.columns
    )
    ovary_summary = summary[
        summary["Study"] == "Ovary (reprocessed), blood-meal time course"
    ]
    assert all(value.endswith("/33") for value in ovary_summary["Samples ≥1 TPM"])

    neuro_plot_before_coloring = _plotly_specs_by_type(app, "box")[1]
    original_condition_order = neuro_plot_before_coloring["layout"]["yaxis"][
        "categoryarray"
    ]
    assert neuro_plot_before_coloring["layout"].get("annotations", []) == []
    point_colors_before = {
        trace["name"]: trace["marker"]["color"]
        for trace in neuro_plot_before_coloring["data"]
        if trace["type"] == "box"
    }
    neuro_color = next(
        widget
        for widget in app.selectbox
        if widget.key == "gene_label_color_neuro_ru"
    )
    neuro_color.set_value("tissue").run()
    assert not app.exception
    neuro_plot = _plotly_specs_by_type(app, "box")[1]
    point_colors_after = {
        trace["name"]: trace["marker"]["color"]
        for trace in neuro_plot["data"]
        if trace["type"] == "box"
    }
    assert point_colors_after == point_colors_before
    scatter_annotations = neuro_plot["layout"]["annotations"]
    colored_condition_order = neuro_plot["layout"]["yaxis"]["categoryarray"]
    assert {annotation["text"] for annotation in scatter_annotations} == set(
        colored_condition_order
    )
    annotation_colors = {
        annotation["text"]: annotation["font"]["color"]
        for annotation in scatter_annotations
    }
    ordered_colors = [annotation_colors[condition] for condition in colored_condition_order]
    color_blocks = [
        color
        for index, color in enumerate(ordered_colors)
        if index == 0 or color != ordered_colors[index - 1]
    ]
    assert len(color_blocks) == len(set(ordered_colors)) == 10
    for color in color_blocks:
        original_positions = [
            original_condition_order.index(condition)
            for condition in colored_condition_order
            if annotation_colors[condition] == color
        ]
        assert original_positions == sorted(original_positions)
    assert neuro_plot["layout"]["yaxis"]["tickfont"]["color"] == "rgba(0,0,0,0)"

    neuro_heatmap = _plotly_specs_by_type(app, "heatmap")[1]
    assert neuro_heatmap["data"][0]["x"] == colored_condition_order
    heatmap_annotations = neuro_heatmap["layout"]["annotations"]
    assert {annotation["text"] for annotation in heatmap_annotations} == set(
        neuro_heatmap["data"][0]["x"]
    )
    assert {
        annotation["text"]: annotation["font"]["color"]
        for annotation in heatmap_annotations
    } == {
        annotation["text"]: annotation["font"]["color"]
        for annotation in scatter_annotations
    }
    assert neuro_heatmap["layout"]["xaxis"]["tickfont"]["color"] == "rgba(0,0,0,0)"

    studies = next(widget for widget in app.multiselect if widget.label == "Studies")
    studies.set_value(["elife", "neuro_ru", "midgut"]).run()
    assert not app.exception
    assert next(
        widget
        for widget in app.selectbox
        if widget.key == "gene_label_color_midgut"
    ).value is None
    assert {
        widget.key
        for widget in app.multiselect
        if str(widget.key).startswith("gene_filter_midgut_")
    } == {
        "gene_filter_midgut_sex",
        "gene_filter_midgut_timepoint",
        "gene_filter_midgut_condition_label",
    }


def test_gene_results_show_aliases_and_missing_studies(monkeypatch):
    monkeypatch.syspath_prepend(str(APP.parent))
    app = AppTest.from_file(str(APP), default_timeout=45).run()
    _select_page(app, "Genes")
    query = next(widget for widget in app.text_input if widget.label == "Genes or identifiers")

    query.set_value("Orco").run()
    assert not app.exception
    matched_title = next(
        element.value
        for element in app.markdown
        if '<div class="matched-genes-title"' in element.value
    )
    assert "Matched genes: 1" in matched_title
    matched_genes = _matched_gene_editor(app).value
    assert matched_genes["Gene"].tolist() == ["Orco"]
    aliases = matched_genes["Alternative names"].item()
    assert "AAEL005776" in aliases
    assert "AaegOr7" in aliases
    assert "Or7" in aliases
    assert "gene14494" in aliases
    single_gene_plots = _plotly_specs_by_type(app, "box")
    assert len(single_gene_plots) == 2
    assert all(_scatter_genes(plot) == {"Orco"} for plot in single_gene_plots)
    single_gene_heatmaps = _plotly_specs_by_type(app, "heatmap")
    assert len(single_gene_heatmaps) == 2

    query.set_value("Ir7a").run()
    assert not app.exception
    query = next(widget for widget in app.text_input if widget.label == "Genes or identifiers")
    assert query.value == "ir7a"
    warnings = " ".join(element.value for element in app.warning)
    assert "Gene not found:" in warnings
    assert "ir7a" in warnings
    assert "Atlas (paper), neurotranscriptome AaegL.RU" in warnings
    matched_genes = _matched_gene_editor(app).value
    assert len(_mean_tpm_columns(matched_genes)) == 2
    assert any(
        matched_genes[column].isna().all()
        for column in _mean_tpm_columns(matched_genes)
    )
    assert matched_genes["Alternative names"].item() == "No alternative names found"

    query.set_value("definitely_not_a_gene").run()
    assert not app.exception
    warnings = [element.value for element in app.warning]
    assert len(warnings) == 2
    assert all("Gene not found:" in warning for warning in warnings)


def test_gene_matched_table_toggle_filters_every_result(monkeypatch):
    monkeypatch.syspath_prepend(str(APP.parent))
    app = AppTest.from_file(str(APP), default_timeout=45).run()
    _select_page(app, "Genes")

    _set_matched_gene_enabled(app, "Ir25a", False)
    assert not app.exception
    summary = next(
        frame.value
        for frame in app.dataframe
        if "Samples ≥1 TPM" in frame.value.columns
    )
    assert set(summary["Gene"]) == {"Orco"}
    raw = next(
        frame.value
        for frame in app.dataframe
        if {"Gene", "Sample", "TPM"} <= set(frame.value.columns)
    )
    assert set(raw["Gene"]) == {"Orco"}
    assert all(
        "Ir25a" not in _scatter_genes(plot)
        for plot in _plotly_specs_by_type(app, "box")
    )
    assert all(
        "Ir25a" not in heatmap["data"][0]["y"]
        for heatmap in _plotly_specs_by_type(app, "heatmap")
    )
    assert any(
        caption.value == "1 of 2 matched genes included in the results below."
        for caption in app.caption
    )


def test_gene_input_normalizes_and_color_codes_submitted_tokens(monkeypatch):
    monkeypatch.syspath_prepend(str(APP.parent))
    app = AppTest.from_file(str(APP), default_timeout=45).run()
    _select_page(app, "Genes")
    query = next(widget for widget in app.text_input if widget.label == "Genes or identifiers")

    query.set_value("oRco, ir7a dlkvnsdlknv").run()
    assert not app.exception
    query = next(widget for widget in app.text_input if widget.label == "Genes or identifiers")
    assert query.value == "orco ir7a dlkvnsdlknv"
    token_html = next(
        element.value
        for element in app.markdown
        if '<div class="gene-token-overlay"' in element.value
    )
    assert "Orco" in token_html
    assert ">OR<" in token_html
    assert "Ir7a" in token_html
    assert ">IR<" in token_html
    assert "dlkvnsdlknv" in token_html
    assert ">not found<" in token_html
    assert token_html.count("gene-token-found") == 2
    assert token_html.count("gene-token-missing") == 1
    assert "gene-token-gene-0" in token_html
    assert "gene-token-gene-6" in token_html
    assert "gene-family-or" in token_html
    assert "gene-family-ir" in token_html
    assert "Click the field to edit as plain text" in token_html
    assert "identifies this gene" in token_html
    assert "Gray means no matching gene was detected" in token_html


def test_family_mode_shows_all_enabled_genes_ranked_by_mean_tpm(monkeypatch):
    monkeypatch.syspath_prepend(str(APP.parent))
    app = AppTest.from_file(str(APP), default_timeout=45).run()
    _select_page(app, "Families")

    assert not app.exception
    source = APP.read_text()
    assert ".st-key-family_setup_panel" in source
    captions = " ".join(element.value for element in app.caption)
    assert "does not combine the family into one score" not in captions
    family_selector = next(
        selectbox for selectbox in app.selectbox if selectbox.label == "Gene family"
    )
    assert family_selector.options[-1] == "Custom family"
    matched_genes = _matched_gene_editor(app).value
    assert len(matched_genes) > 10
    assert matched_genes["Include"].all()
    mean_tpm_columns = _mean_tpm_columns(matched_genes)
    assert len(mean_tpm_columns) == 2
    assert matched_genes[mean_tpm_columns[0]].dropna().is_monotonic_decreasing
    assert all(
        slider.label != "Genes in each heatmap" for slider in app.select_slider
    )
    coverage = next(
        frame.value
        for frame in app.dataframe
        if "Family genes" in frame.value.columns
    )
    assert all(
        widget.label != "Exploratory detection threshold (TPM)"
        for widget in app.number_input
    )
    all_gene_tables = _family_gene_tables(app)
    assert [len(table) for table in all_gene_tables] == coverage[
        "Family genes"
    ].tolist()
    zscore_toggle = next(
        toggle
        for toggle in app.toggle
        if toggle.label == "Show relative pattern within each gene"
    )
    assert "z-scores" in zscore_toggle.help

    heatmaps = [
        _plotly_spec(app, index)
        for index in range(len(app.get("plotly_chart")))
        if _plotly_spec(app, index)["data"][0]["type"] == "heatmap"
    ]
    assert len(all_gene_tables) == len(heatmaps) == 2
    for table, heatmap in zip(all_gene_tables, heatmaps):
        assert table["Mean TPM"].is_monotonic_decreasing
        assert heatmap["data"][0]["y"] == table["Gene"].tolist()


def test_custom_family_reuses_gene_editor_and_shows_all_matches(monkeypatch):
    monkeypatch.syspath_prepend(str(APP.parent))
    app = AppTest.from_file(str(APP), default_timeout=45).run()
    _select_page(app, "Families")

    family_selector = next(
        selectbox for selectbox in app.selectbox if selectbox.label == "Gene family"
    )
    family_selector.set_value("Custom family").run()
    query = next(
        widget for widget in app.text_input if widget.label == "Genes or identifiers"
    )
    query.set_value("Ir25a Orco").run()

    assert not app.exception
    token_html = next(
        element.value
        for element in app.markdown
        if '<div class="gene-token-overlay"' in element.value
    )
    assert "Ir25a" in token_html
    assert "Orco" in token_html
    assert token_html.count("gene-token-found") == 2
    assert all(
        slider.label != "Genes in each heatmap" for slider in app.select_slider
    )
    displayed_tables = _family_gene_tables(app)
    assert displayed_tables
    assert all(len(table) == 2 for table in displayed_tables)
    matched_genes = _matched_gene_editor(app).value
    assert set(matched_genes["Gene"]) == {"Ir25a", "Orco"}
    assert len(_mean_tpm_columns(matched_genes)) == 2

    _set_matched_gene_enabled(app, "Ir25a", False)
    assert not app.exception
    displayed_tables = _family_gene_tables(app)
    assert displayed_tables
    assert all(table["Gene"].tolist() == ["Orco"] for table in displayed_tables)
    source = APP.read_text()
    assert ".st-key-family_query_editor" in source


def test_predefined_family_matched_table_toggle_filters_family(monkeypatch):
    monkeypatch.syspath_prepend(str(APP.parent))
    app = AppTest.from_file(str(APP), default_timeout=45).run()
    _select_page(app, "Families")

    matched_genes = _matched_gene_editor(app).value
    disabled_gene = matched_genes.loc[
        matched_genes["Study coverage"] == "2/2 studies", "Gene"
    ].iloc[0]
    _set_matched_gene_enabled(app, disabled_gene, False)

    assert not app.exception
    displayed_tables = _family_gene_tables(app)
    assert displayed_tables
    assert all(disabled_gene not in table["Gene"].tolist() for table in displayed_tables)


def test_differential_expression_uses_bundled_nfcore_results(monkeypatch):
    monkeypatch.syspath_prepend(str(APP.parent))
    app = AppTest.from_file(str(APP), default_timeout=45).run()
    _select_page(app, "Differential expression")

    study = next(
        selectbox for selectbox in app.selectbox if selectbox.label == "Study"
    )
    study.set_value("midgut").run()

    assert not app.exception, [exception.message for exception in app.exception]
    assert ".st-key-differential_setup_panel" in APP.read_text()
    result_tables = [frame.value for frame in app.dataframe if "FDR" in frame.value.columns]
    assert len(result_tables) == 1
    assert {
        "Gene",
        "DESeq2 base mean",
        "Log₂ fold change",
        "Fold change (target / reference)",
        "Log₂ fold-change SE",
        "Raw p-value",
        "FDR",
        "FDR < 0.05",
    } <= set(result_tables[0].columns)
    assert 10_000 < len(result_tables[0]) <= 19_920

    plot = _plotly_spec(app)
    assert [trace["name"] for trace in plot["data"]] == [
        "FDR ≥ 0.05",
        "FDR < 0.05",
    ]
    assert all(trace["type"] == "scattergl" for trace in plot["data"])
    assert [trace["marker"] for trace in plot["data"]] == [
        {"color": "#66706f", "opacity": 1.0, "size": 3},
        {"color": "#f5b85b", "opacity": 1.0, "size": 3},
    ]
    assert plot["layout"]["xaxis"]["title"]["text"] == "DESeq2 base mean (logarithmic scale)"
    assert plot["layout"]["xaxis"]["type"] == "log"
    assert plot["layout"]["xaxis"]["range"][0] < plot["layout"]["xaxis"]["range"][1]
    assert {"1", "10", "100", "1,000"} <= set(plot["layout"]["xaxis"]["ticktext"])
    assert (
        plot["layout"]["yaxis"]["title"]["text"]
        == "Log₂ fold change (target / reference)"
    )
    assert plot["layout"]["yaxis"]["range"][0] == -plot["layout"]["yaxis"]["range"][1]
    assert plot["layout"]["shapes"][0]["y0"] == 0
    assert plot["layout"]["shapes"][0]["y1"] == 0
    captions = " ".join(element.value for element in app.caption)
    assert "DESeq2" in captions
    assert "ashr-shrunk" in captions
    assert "STAR + Salmon length-scaled gene counts" in captions
    assert "Methods page" in captions
    assert "Significant genes are gold and drawn last" in captions
    assert "Welch" not in captions
    target = next(
        selectbox for selectbox in app.selectbox if selectbox.label == "Target"
    )
    reference = next(
        selectbox for selectbox in app.selectbox if selectbox.label == "Reference"
    )
    assert len(target.options) == 8
    assert len(reference.options) == 7
    assert target.value != reference.value
    original_target = target.value
    original_reference = reference.value
    original_fold_change = result_tables[0]["Log₂ fold change"].to_numpy()

    target.set_value(original_reference).run()
    reference = next(
        selectbox for selectbox in app.selectbox if selectbox.label == "Reference"
    )
    reference.set_value(original_target).run()
    assert not app.exception, [exception.message for exception in app.exception]
    reversed_table = next(
        frame.value for frame in app.dataframe if "FDR" in frame.value.columns
    )
    assert np.allclose(
        reversed_table["Log₂ fold change"], -original_fold_change, equal_nan=True
    )

    fdr_threshold = next(
        widget for widget in app.number_input if widget.label == "FDR threshold"
    )
    assert fdr_threshold.value == 0.05
    fdr_threshold.set_value(0.1).run()
    assert any("FDR < 0.1" in frame.value.columns for frame in app.dataframe)
    threshold_plot = _plotly_spec(app)
    assert [trace["name"] for trace in threshold_plot["data"]] == [
        "FDR ≥ 0.1",
        "FDR < 0.1",
    ]
    assert any(
        button.label == "Download all differential-expression results"
        for button in app.download_button
    )


def test_ovary_differential_expression_has_every_pair(monkeypatch):
    monkeypatch.syspath_prepend(str(APP.parent))
    app = AppTest.from_file(str(APP), default_timeout=45).run()
    _select_page(app, "Differential expression")

    study = next(
        selectbox for selectbox in app.selectbox if selectbox.label == "Study"
    )
    study.set_value("elife").run()

    assert not app.exception, [exception.message for exception in app.exception]
    target = next(
        selectbox for selectbox in app.selectbox if selectbox.label == "Target"
    )
    reference = next(
        selectbox for selectbox in app.selectbox if selectbox.label == "Reference"
    )
    assert len(target.options) == 11
    assert len(reference.options) == 10
    result_tables = [
        frame.value for frame in app.dataframe if "FDR" in frame.value.columns
    ]
    assert len(result_tables) == 1
    assert 10_000 < len(result_tables[0]) <= 19_920


def test_neurotranscriptome_differential_expression_is_not_available(monkeypatch):
    monkeypatch.syspath_prepend(str(APP.parent))
    app = AppTest.from_file(str(APP), default_timeout=45).run()
    _select_page(app, "Differential expression")

    study = next(
        selectbox for selectbox in app.selectbox if selectbox.label == "Study"
    )
    study.set_value("neuro_ru").run()

    assert not app.exception, [exception.message for exception in app.exception]
    rendered_text = " ".join(
        element.value for element in [*app.markdown, *app.info]
    )
    assert "NOT AVAILABLE" in rendered_text
    assert "never computes differential-expression statistics from TPM" in rendered_text
    assert "only for our reprocessed datasets" in rendered_text
    assert not any(
        selectbox.label in {"Target", "Reference"} for selectbox in app.selectbox
    )


def test_cluster_mode_renders_sample_pca(monkeypatch):
    monkeypatch.syspath_prepend(str(APP.parent))
    app = AppTest.from_file(str(APP), default_timeout=45).run()
    _select_page(app, "Clusters")

    assert not app.exception, [exception.message for exception in app.exception]
    assert ".st-key-cluster_setup_panel" in APP.read_text()
    method = _widgets_with_options(app, ["PCA", "UMAP", "t-SNE"])
    assert len(method) == 1
    assert method[0].value == "PCA"
    genes_used = next(
        slider for slider in app.select_slider if slider.label == "Genes used"
    )
    assert genes_used.options[-1] == f"All ({genes_used.value:,})"
    assert "All genes are used by default" in genes_used.help

    plot = _plotly_spec(app)
    assert all(trace["type"] == "scatter" for trace in plot["data"])
    assert plot["layout"]["xaxis"]["title"]["text"].startswith("PC1 (")
    assert plot["layout"]["yaxis"]["title"]["text"].startswith("PC2 (")
    assert plot["layout"]["yaxis"]["showgrid"] is False
    assert plot["layout"]["yaxis"]["zeroline"] is False
    captions = " ".join(element.value for element in app.caption)
    assert "Each point is one biological sample" in captions
    assert "most-variable genes" not in captions
    info_html = next(
        element.value
        for element in app.markdown
        if "Cluster map" in element.value and "section-info-icon" in element.value
    )
    assert f"PCA: {genes_used.value:,} genes (all)" in info_html
    assert "TPM was transformed as log₂(TPM + 1)" in info_html
    assert "PCA preserves broad linear variation" in info_html
    assert "UMAP and t-SNE emphasize local neighborhoods" in info_html

    genes_used.set_value(2_000).run()
    reduced_info_html = next(
        element.value
        for element in app.markdown
        if "Cluster map" in element.value and "section-info-icon" in element.value
    )
    assert "PCA: 2,000 most-variable genes" in reduced_info_html
