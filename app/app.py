from __future__ import annotations

import gzip
import io
import os
from pathlib import Path
import re

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from expression_explorer.comparison import compare_conditions
from expression_explorer.data import (
    DATASET_ORDER,
    expression_long,
    family_members,
    gene_statistics,
    load_datasets,
    matrix_for_genes,
    search_genes,
)


APP_DIR = Path(__file__).resolve().parent
EXPRESSION_DIR = APP_DIR.parent / "expression"
DATA_SCHEMA_VERSION = "2026-07-22-simple-ui-v3"
PUBLIC_MODE = os.getenv("RNA_ATLAS_PUBLIC", "0") == "1"

FAMILIES = {
    "IR · Ionotropic receptors": "Ionotropic receptors (IR)",
    "OR · Odorant receptors": "Odorant receptors (OR)",
    "GR · Gustatory receptors": "Gustatory receptors (GR)",
    "OBP · Odorant-binding proteins": "Odorant-binding proteins (OBP)",
}

ANNOTATION_COLUMNS = [
    "display_name",
    "stable_id",
    "internal_id",
    "raw_symbol",
    "family",
    "drosophila_ortholog",
    "drosophila_blastx_hits",
    "orthodb_category",
    "naming_evidence",
    "search_text",
]

st.set_page_config(
    page_title="Aedes RNA Atlas",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed",
)

st.markdown(
    """
    <style>
    .block-container { padding-top: 2.4rem; padding-bottom: 3rem; max-width: 1320px; }
    .atlas-subtle { color: #9aa8a5; max-width: 820px; margin-bottom: .8rem; }
    .study-note { color: #9aa8a5; font-size: .9rem; }
    [data-testid="stExpander"] { border-color: rgba(128,128,128,.18); }
    </style>
    """,
    unsafe_allow_html=True,
)


@st.cache_resource(show_spinner="Loading expression data…")
def datasets_resource(schema_version: str):
    del schema_version
    return load_datasets(EXPRESSION_DIR)


datasets = datasets_resource(DATA_SCHEMA_VERSION)
ordered_dataset_keys = [key for key in DATASET_ORDER if key in datasets] + sorted(
    key for key in datasets if key not in DATASET_ORDER
)
study_keys = [key for key in ordered_dataset_keys if key != "neuro_legacy"]


def parse_queries(raw: str) -> list[str]:
    values = [item.strip() for item in re.split(r"[,;\n]+", raw) if item.strip()]
    return list(dict.fromkeys(values))


def resolve_one(dataset, query: str) -> pd.DataFrame:
    exact = search_genes(dataset, query, "exact")
    return exact if not exact.empty else search_genes(dataset, query, "contains")


def default_grouping(dataset) -> tuple[str, str]:
    if dataset.key == "elife":
        return "reproductive_state", "Reproductive state"
    if dataset.key.startswith("neuro_"):
        return "tissue_condition", "Tissue + condition"
    return "sample", "Sample"


def grouped_median(long: pd.DataFrame, field: str) -> pd.DataFrame:
    order = [value for value in long[field].dropna().astype(str).unique() if value]
    grouped = long.groupby(["gene", field], as_index=False, sort=False)["tpm"].median()
    grouped[field] = pd.Categorical(grouped[field], categories=order, ordered=True)
    return grouped.sort_values([field, "gene"])


def replicate_figure(
    long: pd.DataFrame,
    field: str,
    field_label: str,
    sort_by_expression: bool = False,
    show_medians: bool = True,
    show_guides: bool = True,
) -> go.Figure:
    plot = long.copy()
    plot["log_tpm"] = np.log2(plot["tpm"] + 1.0)
    plot[field] = plot[field].fillna("Unspecified").astype(str)
    condition_order = plot[field].drop_duplicates().tolist()
    if sort_by_expression:
        condition_order = (
            plot.groupby(field, sort=False)["tpm"]
            .median()
            .sort_values(ascending=False, kind="stable")
            .index.tolist()
        )
    gene_count = plot["gene"].nunique()
    figure = px.strip(
        plot,
        x="log_tpm",
        y=field,
        color="gene" if gene_count > 1 else None,
        orientation="h",
        hover_data={"sample": True, "tpm": ":.3f", "log_tpm": False},
        labels={field: field_label, "log_tpm": "log₂(TPM + 1)", "gene": "Gene"},
        category_orders={field: condition_order},
        color_discrete_sequence=px.colors.qualitative.Set2,
    )
    medians = (
        plot.groupby(["gene", field], as_index=False, sort=False)["log_tpm"]
        .median()
        .rename(columns={"log_tpm": "median"})
    )
    if show_medians and gene_count == 1:
        figure.add_trace(
            go.Scatter(
                x=medians["median"],
                y=medians[field],
                mode="markers",
                name="Group median",
                showlegend=False,
                marker={"symbol": "diamond", "size": 11, "color": "#f5b85b"},
                hovertemplate="%{y}<br>Median log₂(TPM + 1): %{x:.2f}<extra></extra>",
            )
        )
    expression_max = float(plot["log_tpm"].max()) if not plot.empty else 0.0
    expression_upper = max(1.0, expression_max + max(0.2, expression_max * 0.04))
    if show_guides:
        for condition in condition_order[::2]:
            figure.add_shape(
                type="line",
                xref="paper",
                x0=0,
                x1=1,
                yref="y",
                y0=condition,
                y1=condition,
                layer="below",
                line={"color": "rgba(148, 163, 184, 0.28)", "dash": "dot", "width": 1},
            )
    figure.update_traces(jitter=0.34, marker={"opacity": 0.7}, selector={"type": "box"})
    figure.update_layout(
        height=max(340, 110 + 27 * len(condition_order)),
        margin={"l": 20, "r": 20, "t": 20, "b": 45},
        xaxis={
            "title": "log₂(TPM + 1)",
            "range": [0, expression_upper],
            "autorange": False,
        },
        yaxis={
            "title": field_label,
            "categoryorder": "array",
            "categoryarray": condition_order,
            "autorange": "reversed",
            "automargin": True,
        },
        legend={"title": {"text": ""}, "itemclick": False, "itemdoubleclick": False},
    )
    return figure


def heatmap_figure(
    grouped: pd.DataFrame,
    field: str,
    title: str,
    row_zscore: bool,
) -> go.Figure:
    matrix = grouped.pivot(index="gene", columns=field, values="tpm").fillna(0.0)
    transformed = np.log2(matrix + 1.0)
    color_title = "log₂(TPM + 1)"
    colorscale = "Viridis"
    zmid = None
    if row_zscore:
        means = transformed.mean(axis=1)
        stds = transformed.std(axis=1).replace(0, 1.0)
        transformed = transformed.sub(means, axis=0).div(stds, axis=0)
        color_title = "Within-gene z-score"
        colorscale = "RdBu_r"
        zmid = 0
    figure = go.Figure(
        go.Heatmap(
            z=transformed.to_numpy(),
            x=[str(value) for value in transformed.columns],
            y=transformed.index.tolist(),
            colorscale=colorscale,
            zmid=zmid,
            colorbar={"title": color_title},
            hovertemplate="Gene: %{y}<br>Group: %{x}<br>Value: %{z:.2f}<extra></extra>",
        )
    )
    figure.update_layout(
        title=title,
        height=max(390, min(1050, 145 + 23 * len(transformed))),
        margin={"l": 20, "r": 25, "t": 55, "b": 60},
        xaxis={"tickangle": -35},
        yaxis={"autorange": "reversed"},
    )
    return figure


def ma_ratio_log_range(results: pd.DataFrame) -> list[float]:
    finite_ratios = results.loc[
        results["ma_plot_eligible"], "log10_ratio_a_over_b"
    ].dropna()
    if finite_ratios.empty:
        return [-1.0, 1.0]
    limit = max(1.0, float(np.ceil(finite_ratios.abs().quantile(0.995))))
    return [-limit, limit]


def ma_abundance_range(results: pd.DataFrame) -> list[float]:
    finite_abundance = results.loc[
        results["ma_plot_eligible"], "average_tpm"
    ].dropna()
    if finite_abundance.empty:
        return [-1.0, 1.0]
    lower = float(np.floor(np.log10(finite_abundance.quantile(0.005))))
    upper = float(np.ceil(np.log10(finite_abundance.max())))
    return [lower, max(lower + 1.0, upper)]


def base10_ticks(log_range: list[float], suffix: str = "") -> tuple[list[float], list[str]]:
    exponents = range(int(np.ceil(log_range[0])), int(np.floor(log_range[1])) + 1)
    values = [10.0**exponent for exponent in exponents]
    labels = []
    for exponent, value in zip(exponents, values):
        if exponent >= 0:
            label = f"{int(value):,}"
        else:
            label = f"{value:.{-exponent}f}"
        labels.append(f"{label}{suffix}")
    return values, labels


def ma_figure(results: pd.DataFrame, fdr_threshold: float) -> go.Figure:
    plotted = results[results["ma_plot_eligible"]].copy()
    ratio_range = ma_ratio_log_range(results)
    abundance_range = ma_abundance_range(results)
    abundance_ticks, abundance_labels = base10_ticks(abundance_range)
    ratio_ticks, ratio_labels = base10_ticks(ratio_range, "×")
    plotted["passes_fdr"] = plotted["fdr"] < fdr_threshold
    figure = go.Figure()
    for passes_fdr, label, color in (
        (False, f"FDR ≥ {fdr_threshold:g}", "#66706f"),
        (True, f"FDR < {fdr_threshold:g}", "#f5b85b"),
    ):
        subset = plotted[plotted["passes_fdr"].eq(passes_fdr)]
        figure.add_trace(
            go.Scattergl(
                x=subset["average_tpm"],
                y=subset["tpm_ratio_a_over_b"],
                mode="markers",
                name=label,
                marker={"size": 3, "color": color, "opacity": 1.0},
                customdata=subset[
                    [
                        "gene",
                        "stable_id",
                        "mean_tpm_a",
                        "mean_tpm_b",
                        "average_tpm",
                        "fdr",
                    ]
                ].to_numpy(),
                hovertemplate=(
                    "<b>%{customdata[0]}</b><br>"
                    "Stable ID: %{customdata[1]}<br>"
                    "Mean TPM (A): %{customdata[2]:.3f}<br>"
                    "Mean TPM (B): %{customdata[3]:.3f}<br>"
                    "Average TPM: %{customdata[4]:.3f}<br>"
                    "TPM ratio (A / B): %{y:.3f}×<br>"
                    "FDR: %{customdata[5]:.3g}<extra></extra>"
                ),
            )
        )
    figure.add_hline(
        y=1,
        line={"color": "rgba(148,163,184,.38)", "dash": "dot", "width": 1},
    )
    figure.update_layout(
        height=510,
        margin={"l": 25, "r": 25, "t": 35, "b": 55},
        xaxis={
            "title": "Average TPM (logarithmic scale)",
            "type": "log",
            "range": abundance_range,
            "tickmode": "array",
            "tickvals": abundance_ticks,
            "ticktext": abundance_labels,
        },
        yaxis={
            "title": "TPM ratio A / B (logarithmic scale)",
            "type": "log",
            "range": ratio_range,
            "tickmode": "array",
            "tickvals": ratio_ticks,
            "ticktext": ratio_labels,
            "zeroline": False,
        },
        legend={"title": {"text": ""}},
    )
    return figure


def annotation_table(matches: pd.DataFrame) -> pd.DataFrame:
    return matches.reindex(columns=ANNOTATION_COLUMNS, fill_value="").rename(
        columns={
            "display_name": "Gene",
            "stable_id": "Stable ID",
            "internal_id": "Internal ID",
            "raw_symbol": "Published symbol",
            "family": "Family",
            "drosophila_ortholog": "Drosophila ortholog",
            "drosophila_blastx_hits": "Drosophila BLASTX hits",
            "orthodb_category": "OrthoDB category",
            "naming_evidence": "Naming evidence",
            "search_text": "Known aliases",
        }
    )


def download_tsv(label: str, frame: pd.DataFrame, filename: str, key: str):
    st.download_button(
        label,
        frame.to_csv(sep="\t", index=False).encode("utf-8"),
        filename,
        "text/tab-separated-values",
        key=key,
    )


@st.dialog("Import nf-core/rnaseq TPM")
def import_nfcore_dialog():
    st.caption(
        "Use a merged gene-level TPM file such as salmon.merged.gene_tpm.tsv or rsem.merged.gene_tpm.tsv."
    )
    uploaded = st.file_uploader("Merged gene TPM table", type=["tsv", "gz"])
    import_name = st.text_input("Dataset name", value="my_nfcore_run")
    if uploaded is None:
        return
    try:
        payload = uploaded.getvalue()
        compressed = uploaded.name.casefold().endswith(".gz")
        text = gzip.decompress(payload).decode("utf-8") if compressed else payload.decode("utf-8")
        preview = pd.read_csv(io.StringIO(text), sep="\t", nrows=8)
        if preview.shape[1] < 2:
            raise ValueError("Expected one gene-ID column followed by sample columns.")
        numeric = preview.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
        if numeric.isna().all(axis=None):
            raise ValueError("No numeric sample columns were found.")
        st.dataframe(preview, hide_index=True, width="stretch")
        if st.button("Add dataset", type="primary"):
            safe_name = re.sub(r"[^a-zA-Z0-9_-]+", "_", import_name).strip("_") or "nfcore_import"
            suffix = ".merged.gene_tpm.tsv.gz" if compressed else ".merged.gene_tpm.tsv"
            import_dir = EXPRESSION_DIR / "imports"
            import_dir.mkdir(parents=True, exist_ok=True)
            (import_dir / f"{safe_name}{suffix}").write_bytes(payload)
            datasets_resource.clear()
            st.rerun()
    except (OSError, UnicodeError, ValueError, pd.errors.ParserError) as exc:
        st.error(f"Could not read this TPM matrix: {exc}")


with st.sidebar:
    st.markdown("### Data & provenance")
    st.caption("Two published experiments · 155 biological samples")
    if not PUBLIC_MODE and st.button("Import nf-core TPM", width="stretch"):
        import_nfcore_dialog()
    with st.expander("Datasets"):
        for key in study_keys:
            dataset = datasets[key]
            st.markdown(f"**{dataset.label}**")
            st.caption(f"{dataset.paper} · {len(dataset.sample_columns)} sample columns")
        if "neuro_legacy" in datasets:
            st.caption(
                "The AaegL3.3 matrix is retained internally for legacy identifier compatibility; it re-annotates the same 2016 samples and is not a third study."
            )
st.title("Aedes RNA Atlas")
st.markdown(
    '<div class="atlas-subtle">Search genes, compare conditions, or inspect a receptor family. Points are biological samples; diamonds are group medians.</div>',
    unsafe_allow_html=True,
)
st.page_link(
    "pages/1_Mosquito_cheatsheet.py",
    label="Mosquito basics",
    icon="🦟",
)

mode = st.segmented_control(
    "Analysis mode",
    ["Genes", "Families", "Compare conditions"],
    default="Genes",
    label_visibility="collapsed",
)

if mode == "Genes":
    query_text = st.text_input(
        "Genes or identifiers",
        value="Ir25a, Orco",
        help="Comma-separated symbols, AAEL IDs, internal IDs, or aliases.",
    )
    selected_keys = st.multiselect(
        "Studies",
        options=study_keys,
        default=[key for key in ("elife", "neuro_ru") if key in study_keys]
        or study_keys[:1],
        format_func=lambda key: datasets[key].label,
        help="AaegL3.3 is a legacy re-annotation of the same 2016 samples, not a third experiment.",
    )
    sort_column, median_column, guide_column = st.columns(3)
    with sort_column:
        sort_conditions = st.toggle(
            "Sort conditions by expression",
            value=False,
            help="Show the highest median TPM condition first within each study plot.",
        )
    with median_column:
        show_medians = st.toggle(
            "Show group medians",
            value=True,
            help="Show orange diamonds at the group median.",
        )
    with guide_column:
        show_guides = st.toggle(
            "Show row guides",
            value=True,
            help="Show a dotted horizontal guide on alternating condition rows.",
        )
    queries = parse_queries(query_text)

    if not queries:
        st.info("Enter a gene, for example `Ir25a`, `Orco`, or `AAEL005776`.")
    elif not selected_keys:
        st.info("Choose at least one study.")
    else:
        resolved_for_comparison: dict[str, pd.DataFrame] = {}
        unresolved: dict[str, list[str]] = {key: [] for key in selected_keys}

        for query in queries:
            per_study: dict[str, pd.DataFrame] = {}
            for key in selected_keys:
                matches = resolve_one(datasets[key], query)
                if matches.empty:
                    unresolved[key].append(query)
                else:
                    per_study[key] = matches.drop_duplicates("row_id")
                    prior = resolved_for_comparison.get(key)
                    resolved_for_comparison[key] = (
                        per_study[key]
                        if prior is None
                        else pd.concat([prior, per_study[key]], ignore_index=True).drop_duplicates("row_id")
                    )
            if not per_study:
                continue

            first = next(iter(per_study.values())).iloc[0]
            st.markdown(f"## {first['display_name']}")
            subtitle = " · ".join(
                value for value in [first.get("stable_id", ""), first.get("family", "")] if value
            )
            st.caption(subtitle)

            summaries = []
            annotations = []
            raw_rows = []
            for key, matches in per_study.items():
                dataset = datasets[key]
                summary = gene_statistics(dataset, matches)
                summary.insert(0, "Study", dataset.label)
                summary["Samples ≥1 TPM"] = summary.apply(
                    lambda row: f"{round(row['detected_pct'] * len(dataset.sample_columns) / 100):.0f}/{len(dataset.sample_columns)}",
                    axis=1,
                )
                summaries.append(summary)
                annotated = annotation_table(matches)
                annotated.insert(0, "Study", dataset.label)
                annotations.append(annotated)
                long = expression_long(dataset, matches)
                long.insert(0, "Study", dataset.label)
                raw_rows.append(long)

            summary_table = pd.concat(summaries, ignore_index=True).rename(
                columns={
                    "gene": "Gene",
                    "stable_id": "Stable ID",
                    "mean_tpm": "Mean TPM",
                    "median_tpm": "Median TPM",
                    "max_tpm": "Maximum TPM",
                    "top_context": "Top context",
                }
            )
            st.dataframe(
                summary_table[
                    ["Study", "Gene", "Stable ID", "Median TPM", "Maximum TPM", "Samples ≥1 TPM", "Top context"]
                ],
                hide_index=True,
                width="stretch",
                column_config={
                    "Median TPM": st.column_config.NumberColumn(format="%.2f"),
                    "Maximum TPM": st.column_config.NumberColumn(format="%.2f"),
                },
            )

            for key, matches in per_study.items():
                dataset = datasets[key]
                field, field_label = default_grouping(dataset)
                long = expression_long(dataset, matches)
                st.markdown(f"**{dataset.label}**")
                st.plotly_chart(
                    replicate_figure(
                        long,
                        field,
                        field_label,
                        sort_conditions,
                        show_medians,
                        show_guides,
                    ),
                    width="stretch",
                    key=f"gene_plot_{query}_{key}",
                )

            with st.expander("Identifiers, raw values & download"):
                st.dataframe(pd.concat(annotations, ignore_index=True), hide_index=True, width="stretch")
                raw = pd.concat(raw_rows, ignore_index=True)
                raw_table = raw[
                    ["Study", "gene", "stable_id", "sample", "tpm", "tissue", "condition_label", "reproductive_state"]
                ].sort_values(["Study", "gene", "tpm"], ascending=[True, True, False])
                st.dataframe(raw_table, hide_index=True, width="stretch", height=280)
                download_tsv(
                    "Download these values",
                    raw_table,
                    f"{re.sub('[^a-z0-9]+', '_', query.casefold()).strip('_')}_expression.tsv",
                    f"gene_download_{query}",
                )

        if len(queries) > 1 and selected_keys:
            comparisons = []
            for comparison_key in selected_keys:
                comparison_genes = resolved_for_comparison.get(comparison_key)
                if comparison_genes is None or len(comparison_genes) <= 1:
                    continue
                dataset = datasets[comparison_key]
                field, _ = default_grouping(dataset)
                long = expression_long(dataset, comparison_genes)
                grouped = grouped_median(long, field)
                comparisons.append((comparison_key, dataset, field, grouped))

            if comparisons:
                st.markdown("## Compare selected genes")
                st.caption(
                    "One heatmap per study. Colors are scaled within each panel; compare patterns rather than color intensity across papers."
                )
                for comparison_key, dataset, field, grouped in comparisons:
                    st.markdown(f"**{dataset.label}**")
                    st.caption("Group median TPM")
                    st.plotly_chart(
                        heatmap_figure(grouped, field, "", row_zscore=False),
                        width="stretch",
                        key=f"gene_comparison_heatmap_{comparison_key}",
                    )

        missing = {key: values for key, values in unresolved.items() if values}
        if missing:
            with st.expander("Identifiers not found"):
                for key, values in missing.items():
                    st.write(f"**{datasets[key].label}:** {', '.join(values)}")

elif mode == "Families":
    st.caption(
        "Family selection filters to genes annotated as IR, OR, GR, or OBP and ranks individual genes by their highest condition median TPM. It does not combine the family into one score or run a family-level statistical test."
    )
    family_label = st.selectbox("Gene family", list(FAMILIES), index=0)
    family_name = FAMILIES[family_label]
    family_keys = st.multiselect(
        "Studies",
        options=study_keys,
        default=[key for key in ("neuro_ru", "elife") if key in datasets],
        format_func=lambda key: datasets[key].label,
    )
    with st.expander("Advanced display options"):
        top_n = st.slider("Genes in each heatmap", 10, 80, 40, 5)
        threshold = st.number_input("Exploratory detection threshold (TPM)", 0.0, value=1.0, step=0.5)
        row_zscore = st.toggle(
            "Show relative pattern within each gene",
            value=True,
            help="Convert each gene's log₂(TPM + 1) values to z-scores across conditions. Positive means above that gene's average; negative means below it.",
        )

    if not family_keys:
        st.info("Choose at least one study.")
    else:
        coverage = []
        family_data: dict[str, tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, str]] = {}
        for key in family_keys:
            dataset = datasets[key]
            members = family_members(dataset, family_name)
            if members.empty:
                coverage.append({"Study": dataset.label, "Family genes": 0, f"Samples ≥{threshold:g} TPM": "—"})
                continue
            field, _ = default_grouping(dataset)
            long = expression_long(dataset, members)
            group_medians = grouped_median(long, field)
            peak = group_medians.groupby("gene", as_index=False)["tpm"].max().rename(
                columns={"tpm": "Peak group median TPM"}
            )
            detected = (
                long.assign(detected=long["tpm"] >= threshold)
                .groupby("gene", as_index=False)["detected"]
                .sum()
                .rename(columns={"detected": "Detected n"})
            )
            ranking = gene_statistics(dataset, members).merge(peak, on="gene", how="left").merge(
                detected, on="gene", how="left"
            )
            detection_column = f"Samples ≥{threshold:g} TPM"
            ranking[detection_column] = (
                ranking["Detected n"].astype(int).astype(str)
                + f"/{len(dataset.sample_columns)}"
            )
            ranking = ranking.sort_values("Peak group median TPM", ascending=False)

            pinned_names: list[str] = []
            if family_name == "Ionotropic receptors (IR)":
                for pinned in ("Ir25a", "Ir8a", "Ir76b"):
                    if pinned in set(ranking["gene"]):
                        pinned_names.append(pinned)
            ranked_names = [
                gene for gene in ranking["gene"].tolist() if gene not in pinned_names
            ]
            selected_names = ranked_names[: max(0, top_n - len(pinned_names))] + pinned_names
            selected_members = members[members["display_name"].isin(selected_names)]
            selected_long = expression_long(dataset, selected_members)
            selected_grouped = grouped_median(selected_long, field)
            family_data[key] = (members, ranking, selected_grouped, field)
            coverage.append(
                {
                    "Study": dataset.label,
                    "Family genes": len(members),
                    f"Samples ≥{threshold:g} TPM": f"{(ranking['Detected n'] > 0).sum()}/{len(members)} genes",
                }
            )

        st.caption("Family definitions follow each paper's available annotation; studies are not merged.")
        st.dataframe(pd.DataFrame(coverage), hide_index=True, width="stretch")

        for key in family_keys:
            if key not in family_data:
                continue
            dataset = datasets[key]
            members, ranking, selected_grouped, field = family_data[key]
            st.markdown(f"## {dataset.label}")
            st.plotly_chart(
                heatmap_figure(selected_grouped, field, family_label, row_zscore=row_zscore),
                width="stretch",
                key=f"family_heatmap_{key}",
            )
            concise = ranking.rename(
                columns={
                    "gene": "Gene",
                    "stable_id": "Stable ID",
                    "Peak group median TPM": "Peak median TPM",
                    "top_context": "Top context",
                }
            )[
                ["Gene", "Stable ID", "Peak median TPM", detection_column, "Top context"]
            ]
            st.dataframe(
                concise.head(top_n),
                hide_index=True,
                width="stretch",
                column_config={"Peak median TPM": st.column_config.NumberColumn(format="%.2f")},
            )
            with st.expander("Complete family table & download"):
                st.dataframe(concise, hide_index=True, width="stretch", height=340)
                download_tsv(
                    "Download complete family matrix",
                    matrix_for_genes(dataset, members),
                    f"{key}_{family_name.split('(')[-1].rstrip(')').casefold()}_family_tpm.tsv",
                    f"family_download_{key}",
                )

else:
    st.markdown("### Compare all genes between two conditions")
    st.caption(
        "This is an MA plot: right means higher average TPM. A ratio of 1× means equal expression; above 1× means higher in A; below 1× means higher in B."
    )
    comparison_key = st.selectbox(
        "Study",
        options=study_keys,
        format_func=lambda key: datasets[key].label,
    )
    comparison_dataset = datasets[comparison_key]
    comparison_field, _ = default_grouping(comparison_dataset)
    comparison_groups = [
        value
        for value in comparison_dataset.samples[comparison_field]
        .fillna("")
        .astype(str)
        .drop_duplicates()
        if value
    ]
    default_b_name = (
        "6 days post-blood-meal (eggs retained)"
        if comparison_key == "elife"
        else "Antenna · Non-blood-fed / sugar-fed"
    )
    default_b_index = (
        comparison_groups.index(default_b_name)
        if default_b_name in comparison_groups
        else min(1, len(comparison_groups) - 1)
    )
    condition_a_column, condition_b_column = st.columns(2)
    with condition_a_column:
        condition_a = st.selectbox("Condition A", comparison_groups, index=0)
    with condition_b_column:
        condition_b = st.selectbox(
            "Condition B", comparison_groups, index=default_b_index
        )
    fdr_threshold = st.number_input(
        "FDR threshold",
        min_value=0.001,
        max_value=1.0,
        value=0.05,
        step=0.01,
        format="%.3f",
        help="Controls which genes are colored gold, the significant-gene count, and the pass/fail table column.",
    )

    if condition_a == condition_b:
        st.info("Choose two different conditions.")
    else:
        try:
            comparison_results, samples_a, samples_b = compare_conditions(
                comparison_dataset,
                comparison_field,
                condition_a,
                condition_b,
            )
        except ValueError as exc:
            st.warning(str(exc))
        else:
            plotted_results = comparison_results[comparison_results["ma_plot_eligible"]]
            significant_count = int((plotted_results["fdr"] < fdr_threshold).sum())
            omitted_count = len(comparison_results) - len(plotted_results)
            ratio_range = ma_ratio_log_range(comparison_results)
            abundance_range = ma_abundance_range(comparison_results)
            off_scale_count = int(
                (
                    (plotted_results["log10_ratio_a_over_b"] < ratio_range[0])
                    | (plotted_results["log10_ratio_a_over_b"] > ratio_range[1])
                ).sum()
            )
            low_abundance_off_scale = int(
                (np.log10(plotted_results["average_tpm"]) < abundance_range[0]).sum()
            )
            sample_a_metric, sample_b_metric, significant_metric = st.columns(3)
            sample_a_metric.metric("Samples in A", samples_a)
            sample_b_metric.metric("Samples in B", samples_b)
            significant_metric.metric(f"Colored genes · FDR < {fdr_threshold:g}", significant_count)
            st.plotly_chart(
                ma_figure(comparison_results, fdr_threshold),
                width="stretch",
                key=f"condition_comparison_{comparison_key}",
            )
            st.caption(
                f"{len(plotted_results):,} genes plotted. {omitted_count:,} genes with zero mean TPM in A or B are omitted because their A/B ratio is undefined. The initial view excludes {low_abundance_off_scale:,} extreme low-abundance points and {off_scale_count:,} extreme ratios; use Plotly zoom to inspect them."
            )
            st.caption(
                "Gray genes do not pass the selected FDR threshold. Significant genes are gold and drawn last, so gray points cannot cover them. All markers are fully opaque."
            )
            st.caption(
                "Welch's t-test is run on log-transformed replicate TPM; FDR is Benjamini–Hochberg correction across all genes. This is exploratory because TPM-based tests do not model RNA-seq count dispersion. Use raw counts with DESeq2 or edgeR for publication-grade differential expression."
            )

            filter_text = st.text_input(
                "Filter results by gene or Stable ID",
                placeholder="e.g. Ir25a or AAEL005776",
            ).strip()
            displayed_results = comparison_results
            if filter_text:
                needle = filter_text.casefold()
                displayed_results = comparison_results[
                    comparison_results["gene"].astype(str).str.casefold().str.contains(needle, regex=False)
                    | comparison_results["stable_id"].astype(str).str.casefold().str.contains(needle, regex=False)
                ]
            threshold_column = f"FDR < {fdr_threshold:g}"
            displayed_results = displayed_results.assign(
                passes_fdr=displayed_results["fdr"] < fdr_threshold
            )
            display_table = displayed_results.rename(
                columns={
                    "gene": "Gene",
                    "stable_id": "Stable ID",
                    "mean_tpm_a": "Mean TPM (A)",
                    "mean_tpm_b": "Mean TPM (B)",
                    "average_tpm": "Average TPM",
                    "median_tpm_a": "Median TPM (A)",
                    "median_tpm_b": "Median TPM (B)",
                    "tpm_ratio_a_over_b": "TPM ratio (A / B)",
                    "p_value": "Raw p-value",
                    "fdr": "FDR",
                    "passes_fdr": threshold_column,
                }
            )[
                [
                    "Gene",
                    "Stable ID",
                    "Mean TPM (A)",
                    "Mean TPM (B)",
                    "Average TPM",
                    "Median TPM (A)",
                    "Median TPM (B)",
                    "TPM ratio (A / B)",
                    "Raw p-value",
                    "FDR",
                    threshold_column,
                ]
            ]
            st.dataframe(
                display_table,
                hide_index=True,
                width="stretch",
                height=620,
                column_config={
                    "Mean TPM (A)": st.column_config.NumberColumn(format="%.3f"),
                    "Mean TPM (B)": st.column_config.NumberColumn(format="%.3f"),
                    "Average TPM": st.column_config.NumberColumn(format="%.3f"),
                    "Median TPM (A)": st.column_config.NumberColumn(format="%.3f"),
                    "Median TPM (B)": st.column_config.NumberColumn(format="%.3f"),
                    "TPM ratio (A / B)": st.column_config.NumberColumn(format="%.3f"),
                    "Raw p-value": st.column_config.NumberColumn(format="%.3e"),
                    "FDR": st.column_config.NumberColumn(format="%.3e"),
                },
            )

            download_results = comparison_results.assign(
                passes_fdr=comparison_results["fdr"] < fdr_threshold
            )
            download_table = download_results.rename(
                columns={
                    "gene": "Gene",
                    "stable_id": "Stable ID",
                    "mean_tpm_a": "Mean TPM (A)",
                    "mean_tpm_b": "Mean TPM (B)",
                    "average_tpm": "Average TPM",
                    "median_tpm_a": "Median TPM (A)",
                    "median_tpm_b": "Median TPM (B)",
                    "tpm_ratio_a_over_b": "TPM ratio (A / B)",
                    "p_value": "Raw p-value",
                    "fdr": "FDR",
                    "passes_fdr": threshold_column,
                }
            )[
                [
                    "Gene",
                    "Stable ID",
                    "Mean TPM (A)",
                    "Mean TPM (B)",
                    "Average TPM",
                    "Median TPM (A)",
                    "Median TPM (B)",
                    "TPM ratio (A / B)",
                    "Raw p-value",
                    "FDR",
                    threshold_column,
                ]
            ]
            download_table.insert(0, "Condition B", condition_b)
            download_table.insert(0, "Condition A", condition_a)
            download_table.insert(0, "Study", comparison_dataset.label)
            download_tsv(
                "Download all comparison results",
                download_table,
                f"{comparison_key}_condition_comparison.tsv",
                f"condition_comparison_download_{comparison_key}",
            )
