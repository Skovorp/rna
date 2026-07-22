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
    with st.expander("Methods & limits"):
        st.markdown(
            "TPM supports descriptive expression patterns within a study. Studies are plotted separately; raw TPM should not be ranked across papers. Differential-expression claims require raw counts and a replicate-aware model."
        )

st.title("Aedes RNA Atlas")
st.markdown(
    '<div class="atlas-subtle">Search genes, compare a panel, or inspect a receptor family. Points are biological samples; diamonds are group medians.</div>',
    unsafe_allow_html=True,
)
st.page_link(
    "pages/1_Mosquito_cheatsheet.py",
    label="Mosquito anatomy & stages cheatsheet",
    icon="🦟",
)

mode = st.segmented_control(
    "Analysis mode",
    ["Genes", "Families"],
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
        default=["neuro_ru"] if "neuro_ru" in datasets else study_keys[:1],
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
                summary["Detected samples"] = summary.apply(
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
                    ["Study", "Gene", "Stable ID", "Median TPM", "Maximum TPM", "Detected samples", "Top context"]
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
            comparison_key = selected_keys[0]
            comparison_genes = resolved_for_comparison.get(comparison_key)
            if comparison_genes is not None and len(comparison_genes) > 1:
                dataset = datasets[comparison_key]
                field, _ = default_grouping(dataset)
                long = expression_long(dataset, comparison_genes)
                grouped = grouped_median(long, field)
                st.markdown("## Compare selected genes")
                st.caption(f"Group medians · {dataset.label}")
                st.plotly_chart(
                    heatmap_figure(grouped, field, "", row_zscore=False),
                    width="stretch",
                    key="gene_comparison_heatmap",
                )

        missing = {key: values for key, values in unresolved.items() if values}
        if missing:
            with st.expander("Identifiers not found"):
                for key, values in missing.items():
                    st.write(f"**{datasets[key].label}:** {', '.join(values)}")

else:
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
        row_zscore = st.toggle("Show relative pattern within each gene", value=True)

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
            ranking["Detected samples"] = ranking["Detected n"].astype(int).astype(str) + f"/{len(dataset.sample_columns)}"
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
                ["Gene", "Stable ID", "Peak median TPM", "Detected samples", "Top context"]
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

with st.expander("How to interpret these plots"):
    st.markdown(
        "Points are individual samples and diamonds are group medians. TPM is descriptive normalized abundance. Compare patterns within a study; across papers, compare qualitative patterns only. Family z-scores emphasize where each gene is relatively enriched, while the adjacent table preserves absolute TPM context."
    )
