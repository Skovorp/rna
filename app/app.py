from __future__ import annotations

import html
from pathlib import Path
import re

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
import streamlit.components.v1 as components

from expression_explorer.clustering import METHODS, sample_embedding
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
DATA_SCHEMA_VERSION = "2026-08-04-midgut-v1"

FAMILIES = {
    "IR · Ionotropic receptors": "Ionotropic receptors (IR)",
    "OR · Odorant receptors": "Odorant receptors (OR)",
    "GR · Gustatory receptors": "Gustatory receptors (GR)",
    "OBP · Odorant-binding proteins": "Odorant-binding proteins (OBP)",
}
CUSTOM_FAMILY_LABEL = "Custom family"

GENE_X_SCALES = ("Linear", "Log base 2", "Log base 10")
GENE_LOG_BASES = {
    "Log base 2": 2.0,
    "Log base 10": 10.0,
}

GENE_FAMILY_SHORT_LABELS = {
    "Ionotropic receptors (IR)": "IR",
    "Odorant receptors (OR)": "OR",
    "Gustatory receptors (GR)": "GR",
    "Odorant-binding proteins (OBP)": "OBP",
}

GENE_COLOR_NAMES = (
    "blue",
    "amber",
    "pink",
    "cyan",
    "violet",
    "orange",
    "lime",
    "rose",
)

GENE_COLOR_VALUES = (
    "#3b82f6",
    "#f59e0b",
    "#ec4899",
    "#06b6d4",
    "#8b5cf6",
    "#f97316",
    "#84cc16",
    "#f43f5e",
)

GENE_FAMILY_COLOR_NAMES = {
    "IR": "purple",
    "OR": "emerald",
    "GR": "orange",
    "OBP": "pink",
    "Other": "slate",
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
    [data-testid="stHeader"], [data-testid="stToolbar"],
    [data-testid="stSidebar"], [data-testid="stSidebarCollapsedControl"],
    [data-testid="collapsedControl"], [data-testid="stStatusWidget"],
    [data-testid="stAppDeployButton"], #MainMenu, footer {
        display: none !important;
    }
    .block-container { padding-top: 0; padding-bottom: 3rem; max-width: 1320px; }
    .atlas-subtle { color: #9aa8a5; max-width: 820px; margin-bottom: .8rem; }
    .study-note { color: #9aa8a5; font-size: .9rem; }
    .st-key-site_nav {
        position: sticky;
        top: 0;
        z-index: 999;
        margin-bottom: 1.4rem;
        padding: .75rem 0;
        background: rgba(14, 17, 23, .96);
        backdrop-filter: blur(14px);
    }
    .st-key-home_logo button {
        justify-content: flex-start;
        min-height: 2.4rem;
        padding: 0;
        border: 0;
        background: transparent;
        color: inherit;
        font-size: 1.05rem;
        font-weight: 800;
        letter-spacing: -.015em;
        white-space: nowrap;
    }
    .st-key-home_logo button:hover,
    .st-key-home_logo button:focus { color: #6ee7b7; background: transparent; }
    .page-heading { margin: .6rem 0 1.4rem; }
    .page-heading h1 { margin: 0 0 .35rem; font-size: 2.2rem; line-height: 1.08; }
    .page-heading p { margin: 0; max-width: 54rem; color: #9aa8a5; font-size: 1rem; }
    .gene-input-label { font-size: .875rem; font-weight: 600; margin-bottom: .2rem; }
    .st-key-gene_query_editor, .st-key-family_query_editor { position: relative; }
    .st-key-gene_query_editor:not(.gene-query-editing) input,
    .st-key-family_query_editor:not(.gene-query-editing) input {
        color: transparent !important;
        caret-color: transparent !important;
    }
    .st-key-gene_query_editor .gene-token-overlay,
    .st-key-family_query_editor .gene-token-overlay {
        position: absolute;
        z-index: 3;
        top: .28rem;
        left: .7rem;
        right: .7rem;
        height: 2rem;
        display: flex;
        align-items: center;
        gap: .38rem;
        overflow: hidden;
        white-space: nowrap;
        pointer-events: none;
        transition: opacity .08s ease;
    }
    .st-key-gene_query_editor.gene-query-editing .gene-token-overlay,
    .st-key-family_query_editor.gene-query-editing .gene-token-overlay {
        opacity: 0;
        visibility: hidden;
    }
    .st-key-gene_query_editor .stElementContainer:has(.gene-token-overlay),
    .st-key-family_query_editor .stElementContainer:has(.gene-token-overlay) {
        position: absolute !important;
        inset: 0 !important;
        height: 100%;
        min-height: 0;
        margin: 0;
        pointer-events: none;
    }
    .gene-token {
        display: inline-flex;
        align-items: center;
        gap: .3rem;
        flex: 0 0 auto;
        border: 1px solid;
        border-radius: 999px;
        padding: .18rem .55rem;
        color: #17202a;
        font-size: .78rem;
        font-weight: 650;
        line-height: 1.25rem;
    }
    .gene-token-gene-0 { background: #bfdbfe; border-color: #3b82f6; }
    .gene-token-gene-1 { background: #fde68a; border-color: #f59e0b; }
    .gene-token-gene-2 { background: #fbcfe8; border-color: #ec4899; }
    .gene-token-gene-3 { background: #a5f3fc; border-color: #06b6d4; }
    .gene-token-gene-4 { background: #ddd6fe; border-color: #8b5cf6; }
    .gene-token-gene-5 { background: #fed7aa; border-color: #f97316; }
    .gene-token-gene-6 { background: #d9f99d; border-color: #84cc16; }
    .gene-token-gene-7 { background: #fecdd3; border-color: #f43f5e; }
    .gene-token-missing { background: #e5e7eb; border-color: #9ca3af; }
    .gene-token-family {
        border-radius: 999px;
        padding: .03rem .3rem;
        color: white;
        font-size: .65rem;
        font-weight: 800;
        letter-spacing: .015em;
    }
    .gene-family-ir { background: #7c3aed; }
    .gene-family-or { background: #059669; }
    .gene-family-gr { background: #ea580c; }
    .gene-family-obp { background: #db2777; }
    .gene-family-other, .gene-family-mixed, .gene-family-missing {
        background: #64748b;
    }
    .st-key-gene_setup_panel,
    .st-key-family_setup_panel,
    .st-key-condition_setup_panel,
    .st-key-cluster_setup_panel {
        margin: .7rem 0 1rem;
        padding: .9rem 1rem .7rem;
        border: 1px solid rgba(148, 163, 184, .25);
        border-radius: .75rem;
        background: rgba(148, 163, 184, .09);
    }
    .matched-genes-panel {
        margin: .6rem 0 1.2rem;
        border: 1px solid rgba(148, 163, 184, .22);
        border-radius: .75rem;
        overflow: hidden;
    }
    .matched-genes-title {
        padding: .55rem .75rem;
        border-bottom: 1px solid rgba(148, 163, 184, .18);
        background: rgba(148, 163, 184, .08);
        font-size: .82rem;
        font-weight: 750;
    }
    .matched-genes-scroll {
        max-height: 15rem;
        overflow-y: auto;
        overscroll-behavior: contain;
    }
    .matched-gene-row {
        display: grid;
        grid-template-columns: minmax(8rem, 1fr) minmax(14rem, 3fr) auto;
        gap: .75rem;
        align-items: center;
        padding: .5rem .7rem;
        border-bottom: 1px solid rgba(148, 163, 184, .12);
        border-left: 4px solid;
        font-size: .78rem;
    }
    .matched-gene-row:last-child { border-bottom: 0; }
    .matched-gene-name { display: flex; align-items: center; gap: .4rem; font-weight: 750; }
    .matched-gene-aliases { color: #9aa8a5; overflow-wrap: anywhere; }
    .matched-gene-coverage { color: #9aa8a5; white-space: nowrap; }
    .section-title-with-info {
        display: flex;
        align-items: center;
        gap: .5rem;
        margin: 1.8rem 0 .8rem;
    }
    .section-title-with-info h2 { margin: 0; padding: 0; }
    .section-info-icon {
        position: relative;
        display: inline-flex;
        align-items: center;
        justify-content: center;
        width: 1.15rem;
        height: 1.15rem;
        border: 1px solid #9aa8a5;
        border-radius: 999px;
        color: #9aa8a5;
        font-size: .72rem;
        font-weight: 800;
        cursor: help;
    }
    .section-info-tooltip {
        position: absolute;
        z-index: 1000;
        top: 50%;
        left: 1.55rem;
        width: min(28rem, 70vw);
        padding: .55rem .65rem;
        border: 1px solid rgba(148, 163, 184, .35);
        border-radius: .45rem;
        background: #1b2229;
        box-shadow: 0 .45rem 1.4rem rgba(0, 0, 0, .35);
        color: #e5e7eb;
        font-size: .78rem;
        font-weight: 400;
        line-height: 1.4;
        text-align: left;
        opacity: 0;
        visibility: hidden;
        pointer-events: none;
        transform: translateY(-50%);
    }
    .section-info-icon:hover .section-info-tooltip,
    .section-info-icon:focus .section-info-tooltip {
        opacity: 1;
        visibility: visible;
    }
    @media (max-width: 800px) {
        .st-key-site_nav [data-testid="stHorizontalBlock"] { gap: .35rem; }
        .matched-gene-row { grid-template-columns: 1fr; gap: .25rem; }
    }
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


@st.cache_data(show_spinner="Computing sample map…")
def cluster_embedding_cached(
    schema_version: str,
    dataset_key: str,
    method: str,
    variable_genes: int,
):
    dataset = datasets_resource(schema_version)[dataset_key]
    return sample_embedding(dataset, method, variable_genes)


def parse_queries(raw: str) -> list[str]:
    values = [item.strip() for item in re.split(r"[,;\s]+", raw) if item.strip()]
    return list(dict.fromkeys(values))


def normalize_gene_query_input(state_key: str = "gene_query_text") -> None:
    raw = st.session_state.get(state_key, "")
    st.session_state[state_key] = " ".join(
        query.casefold() for query in parse_queries(raw)
    )


def gene_color_index(gene_name: str) -> int:
    return sum(
        (position + 1) * ord(character)
        for position, character in enumerate(gene_name.casefold())
    ) % len(GENE_COLOR_NAMES)


def resolve_one(dataset, query: str) -> pd.DataFrame:
    exact = search_genes(dataset, query, "exact")
    return exact if not exact.empty else search_genes(dataset, query, "contains")


def detected_gene_tokens(
    queries: list[str],
    selected_keys: list[str],
) -> list[dict[str, str]]:
    tokens = []
    for query in queries:
        matched_frames = []
        found_keys = []
        for key in selected_keys:
            matches = resolve_one(datasets[key], query)
            if not matches.empty:
                matched_frames.append(matches)
                found_keys.append(key)
        if not matched_frames:
            tokens.append(
                {
                    "label": query,
                    "family": "not found",
                    "status": "missing",
                    "title": "Not found in the selected studies",
                }
            )
            continue

        matches = pd.concat(matched_frames, ignore_index=True)
        display_names = list(
            dict.fromkeys(
                str(value).strip()
                for value in matches["display_name"]
                if str(value).strip()
            )
        )
        families = list(
            dict.fromkeys(
                GENE_FAMILY_SHORT_LABELS.get(str(value), "Other")
                for value in matches["family"]
            )
        )
        if len(display_names) == 1:
            label = display_names[0]
        else:
            label = f"{query} ({len(display_names)} genes)"
        found_labels = ", ".join(datasets[key].label for key in found_keys)
        tokens.append(
            {
                "label": label,
                "family": "/".join(families),
                "status": "found",
                "title": f"Input: {query} · Found in: {found_labels}",
            }
        )
    return tokens


def render_gene_tokens(tokens: list[dict[str, str]]) -> None:
    if not tokens:
        return
    chips = []
    for token in tokens:
        label = html.escape(token["label"])
        family = html.escape(token["family"])
        status = token["status"]
        if status == "found":
            color_index = gene_color_index(token["label"])
            gene_color = GENE_COLOR_NAMES[color_index]
            family_key = token["family"] if token["family"] in GENE_FAMILY_COLOR_NAMES else "mixed"
            family_color = GENE_FAMILY_COLOR_NAMES.get(family_key, "slate")
            chip_classes = f"gene-token-found gene-token-gene-{color_index}"
            family_class = f"gene-family-{family_key.casefold()}"
            explanation = (
                f"{token['title']}. {gene_color.capitalize()} identifies this gene; "
                f"the {family_color} family tag identifies {token['family']}. "
                "Click the field to edit as plain text, then press Enter to redraw."
            )
        else:
            chip_classes = "gene-token-missing"
            family_class = "gene-family-missing"
            explanation = (
                f"{token['title']}. Gray means no matching gene was detected. "
                "Click the field to edit as plain text, then press Enter to redraw."
            )
        title = html.escape(explanation, quote=True)
        chips.append(
            f'<span class="gene-token {chip_classes}" title="{title}">'
            f"<span>{label}</span>"
            f'<span class="gene-token-family {family_class}">{family}</span>'
            "</span>"
        )
    st.markdown(
        f'<div class="gene-token-overlay" aria-label="Parsed genes">{"".join(chips)}</div>',
        unsafe_allow_html=True,
    )


def matched_gene_entries(
    per_study: dict[str, pd.DataFrame],
) -> list[dict[str, object]]:
    entries: dict[str, dict[str, object]] = {}
    for study_key, matches in per_study.items():
        for gene in matches.itertuples(index=False):
            display_name = str(gene.display_name).strip()
            if not display_name:
                continue
            entry_key = display_name.casefold()
            entry = entries.setdefault(
                entry_key,
                {
                    "display_name": display_name,
                    "families": {},
                    "aliases": {},
                    "study_keys": [],
                },
            )
            family = GENE_FAMILY_SHORT_LABELS.get(str(gene.family), "Other")
            entry["families"].setdefault(family.casefold(), family)
            for alias in gene.aliases:
                alias_text = str(alias).strip()
                if alias_text and alias_text.casefold() != entry_key:
                    entry["aliases"].setdefault(alias_text.casefold(), alias_text)
            if study_key not in entry["study_keys"]:
                entry["study_keys"].append(study_key)
    return list(entries.values())


def render_matched_gene_list(
    entries: list[dict[str, object]],
    selected_study_count: int,
) -> None:
    if not entries:
        return
    rows = []
    for entry in entries:
        display_name = str(entry["display_name"])
        families = list(entry["families"].values())
        family = "/".join(families)
        family_key = family if family in GENE_FAMILY_COLOR_NAMES else "mixed"
        aliases = sorted(entry["aliases"].values(), key=str.casefold)
        alias_text = ", ".join(aliases) if aliases else "No alternative names found"
        found_count = len(entry["study_keys"])
        color = GENE_COLOR_VALUES[gene_color_index(display_name)]
        rows.append(
            f'<div class="matched-gene-row" style="border-left-color:{color}">'
            '<div class="matched-gene-name">'
            f"<span>{html.escape(display_name)}</span>"
            f'<span class="gene-token-family gene-family-{family_key.casefold()}">'
            f"{html.escape(family)}</span></div>"
            f'<div class="matched-gene-aliases">{html.escape(alias_text)}</div>'
            f'<div class="matched-gene-coverage">{found_count}/{selected_study_count} studies</div>'
            "</div>"
        )
    st.markdown(
        '<div class="matched-genes-panel">'
        f'<div class="matched-genes-title">Matched genes · {len(entries)}</div>'
        f'<div class="matched-genes-scroll">{"".join(rows)}</div>'
        "</div>",
        unsafe_allow_html=True,
    )


def install_gene_editor_behavior(editor_key: str, input_key: str) -> None:
    script = """
        <script>
        (() => {
          const doc = window.parent.document;
          const bind = () => {
            const editor = doc.querySelector('__EDITOR_SELECTOR__');
            const input = doc.querySelector('__INPUT_SELECTOR__');
            if (!editor || !input) return false;
            if (input.dataset.geneEditorBound === '__BINDING_KEY__') return true;
            input.dataset.geneEditorBound = '__BINDING_KEY__';

            if (doc.activeElement === input) input.blur();
            editor.classList.remove('gene-query-editing');

            input.addEventListener('focus', () => {
              editor.classList.add('gene-query-editing');
            });
            input.addEventListener('blur', () => {
              editor.classList.remove('gene-query-editing');
            });
            input.addEventListener('keydown', (event) => {
              if (event.key !== 'Enter') return;
              window.setTimeout(() => {
                input.blur();
                editor.classList.remove('gene-query-editing');
              }, 0);
            });
            return true;
          };

          if (!bind()) {
            const observer = new MutationObserver(() => {
              if (bind()) observer.disconnect();
            });
            observer.observe(doc.body, { childList: true, subtree: true });
            window.setTimeout(() => observer.disconnect(), 10000);
          }
        })();
        </script>
        """
    script = (
        script.replace("__EDITOR_SELECTOR__", f".st-key-{editor_key}")
        .replace("__INPUT_SELECTOR__", f".st-key-{input_key} input")
        .replace("__BINDING_KEY__", editor_key)
    )
    components.html(
        script,
        height=0,
        scrolling=False,
    )


def render_gene_query_editor(
    *,
    state_key: str,
    editor_key: str,
    selected_keys: list[str],
    default_value: str = "",
) -> list[str]:
    if state_key not in st.session_state:
        st.session_state[state_key] = default_value
    st.markdown(
        '<div class="gene-input-label">Genes or identifiers</div>',
        unsafe_allow_html=True,
    )
    with st.container(key=editor_key):
        query_text = st.text_input(
            "Genes or identifiers",
            key=state_key,
            on_change=normalize_gene_query_input,
            args=(state_key,),
            label_visibility="collapsed",
            help="Separate symbols or IDs with spaces, commas, semicolons, or new lines. Press Enter to apply.",
        )
        queries = parse_queries(query_text)
        render_gene_tokens(detected_gene_tokens(queries, selected_keys))
    install_gene_editor_behavior(editor_key, state_key)
    return queries


def resolve_queries(dataset, queries: list[str]) -> pd.DataFrame:
    matches = [resolve_one(dataset, query) for query in queries]
    matches = [frame for frame in matches if not frame.empty]
    if not matches:
        return dataset.genes.iloc[0:0].copy()
    return pd.concat(matches, ignore_index=True).drop_duplicates("row_id")


def gene_count_options(maximum: int) -> list[int]:
    """Return slider stops whose final option always includes every gene."""
    if maximum <= 0:
        return [0]
    if maximum < 10:
        return list(range(1, maximum + 1))
    return list(range(10, maximum, 5)) + [maximum]


def alternative_names_by_gene(
    per_study: dict[str, pd.DataFrame],
) -> list[tuple[str, list[str]]]:
    """Collect every non-primary alias found across the selected studies."""
    display_names: dict[str, str] = {}
    alternatives: dict[str, dict[str, str]] = {}
    for matches in per_study.values():
        for _, gene in matches.iterrows():
            display_name = str(gene.get("display_name", "")).strip()
            if not display_name:
                continue
            gene_key = display_name.casefold()
            display_names.setdefault(gene_key, display_name)
            alternatives.setdefault(gene_key, {})
            for alias in gene.get("aliases", ()):
                alias_text = str(alias).strip()
                alias_key = alias_text.casefold()
                if alias_text and alias_key != gene_key:
                    alternatives[gene_key].setdefault(alias_key, alias_text)
    return [
        (
            display_names[gene_key],
            sorted(alternatives[gene_key].values(), key=str.casefold),
        )
        for gene_key in display_names
    ]


def default_grouping(dataset) -> tuple[str, str]:
    if dataset.key == "elife":
        return "reproductive_state", "Reproductive state"
    if dataset.key.startswith("neuro_"):
        return "tissue_condition", "Tissue + condition"
    if dataset.key == "midgut":
        return "condition_label", "Sex + blood-meal time"
    return "sample", "Sample"


def grouped_median(long: pd.DataFrame, field: str) -> pd.DataFrame:
    order = [value for value in long[field].dropna().astype(str).unique() if value]
    grouped = long.groupby(["gene", field], as_index=False, sort=False)["tpm"].median()
    grouped[field] = pd.Categorical(grouped[field], categories=order, ordered=True)
    return grouped.sort_values([field, "gene"])


def tpm_axis_position(
    values: pd.Series | pd.DataFrame, scale: str
) -> pd.Series | pd.DataFrame:
    numeric = values.astype(float)
    if scale == "Linear":
        return numeric
    if scale not in GENE_LOG_BASES:
        raise ValueError(f"Unknown gene x-axis scale: {scale}")
    return np.log1p(numeric) / np.log(GENE_LOG_BASES[scale])


def tpm_axis_config(max_tpm: float, scale: str) -> dict:
    if scale == "Linear":
        upper = max(1.0, max_tpm + max(0.2, max_tpm * 0.04))
        return {
            "title": "TPM (linear scale)",
            "range": [0, upper],
            "autorange": False,
        }
    if scale not in GENE_LOG_BASES:
        raise ValueError(f"Unknown gene x-axis scale: {scale}")

    base = GENE_LOG_BASES[scale]
    max_power = max(0, int(np.ceil(np.log(max(max_tpm, 1.0)) / np.log(base))))
    actual_ticks = [0.0] + [base**power for power in range(max_power + 1)]
    tick_positions = [float(np.log1p(value) / np.log(base)) for value in actual_ticks]
    transformed_max = float(np.log1p(max_tpm) / np.log(base))
    upper = max(transformed_max, tick_positions[-1])
    upper += max(0.08, upper * 0.025)
    base_label = "log₂" if base == 2 else "log₁₀"
    return {
        "title": f"TPM ({base_label} scale)",
        "range": [0, upper],
        "autorange": False,
        "tickmode": "array",
        "tickvals": tick_positions,
        "ticktext": [f"{int(value):,}" for value in actual_ticks],
    }


def replicate_figure(
    long: pd.DataFrame,
    field: str,
    field_label: str,
    sort_by_expression: bool = False,
    show_medians: bool = True,
    show_guides: bool = True,
    x_scale: str = "Log base 2",
) -> go.Figure:
    plot = long.copy()
    plot["axis_tpm"] = tpm_axis_position(plot["tpm"], x_scale)
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
    gene_names = plot["gene"].drop_duplicates().astype(str).tolist()
    gene_color_map = {
        gene_name: GENE_COLOR_VALUES[gene_color_index(gene_name)]
        for gene_name in gene_names
    }
    figure = px.strip(
        plot,
        x="axis_tpm",
        y=field,
        color="gene",
        orientation="h",
        hover_data={"sample": True, "tpm": ":.3f", "axis_tpm": False},
        labels={
            field: field_label,
            "axis_tpm": "TPM",
            "tpm": "TPM",
            "sample": "Sample",
            "gene": "Gene",
        },
        category_orders={field: condition_order},
        color_discrete_map=gene_color_map,
    )
    if gene_count == 1:
        figure.update_traces(marker={"color": gene_color_map[gene_names[0]]})
    medians = (
        plot.groupby(["gene", field], as_index=False, sort=False)["tpm"]
        .median()
    )
    medians["axis_tpm"] = tpm_axis_position(medians["tpm"], x_scale)
    if show_medians:
        trace_colors = {
            str(trace.name): trace.marker.color
            for trace in figure.data
            if getattr(trace, "marker", None) is not None
        }
        if gene_count == 1:
            trace_colors[gene_names[0]] = gene_color_map[gene_names[0]]
        for gene_name, gene_medians in medians.groupby("gene", sort=False):
            figure.add_trace(
                go.Scatter(
                    x=gene_medians["axis_tpm"],
                    y=gene_medians[field],
                    mode="markers",
                    name="Group median",
                    legendgroup=str(gene_name),
                    showlegend=False,
                    marker={
                        "symbol": "diamond",
                        "size": 11,
                        "color": trace_colors.get(str(gene_name), "#f5b85b"),
                        "line": {"color": "#111827", "width": 0.7},
                    },
                    customdata=gene_medians[["gene", "tpm"]].to_numpy(),
                    hovertemplate=(
                        "<b>%{customdata[0]}</b><br>%{y}"
                        "<br>Median TPM: %{customdata[1]:.3f}<extra></extra>"
                    ),
                )
            )
    expression_max = float(plot["tpm"].max()) if not plot.empty else 0.0
    if show_guides:
        for condition in condition_order:
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
        xaxis=tpm_axis_config(expression_max, x_scale),
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
    value_scale: str = "Log base 2",
) -> go.Figure:
    matrix = grouped.pivot(index="gene", columns=field, values="tpm").fillna(0.0)
    transformed = tpm_axis_position(matrix, value_scale)
    color_title = f"Median TPM<br>({value_scale.casefold()} scale)"
    colorscale = "Viridis"
    zmid = None
    zmin = 0
    colorbar: dict = {"title": color_title}
    hovertemplate = (
        "Gene: %{y}<br>Group: %{x}<br>Median TPM: %{customdata:.3f}<extra></extra>"
    )
    if row_zscore:
        transformed = np.log2(matrix + 1.0)
        means = transformed.mean(axis=1)
        stds = transformed.std(axis=1).replace(0, 1.0)
        transformed = transformed.sub(means, axis=0).div(stds, axis=0)
        color_title = "Within-gene z-score"
        colorscale = "RdBu_r"
        zmid = 0
        zmin = None
        colorbar = {"title": color_title}
        hovertemplate = (
            "Gene: %{y}<br>Group: %{x}<br>Within-gene z-score: %{z:.2f}"
            "<br>Median TPM: %{customdata:.3f}<extra></extra>"
        )
    elif value_scale in GENE_LOG_BASES:
        max_tpm = float(matrix.to_numpy().max()) if matrix.size else 0.0
        axis_config = tpm_axis_config(max_tpm, value_scale)
        transformed_max = float(np.asarray(transformed).max()) if matrix.size else 0.0
        visible_ticks = [
            (position, label)
            for position, label in zip(
                axis_config["tickvals"], axis_config["ticktext"]
            )
            if position <= transformed_max + 1e-12
        ]
        colorbar = {
            "title": color_title,
            "tickmode": "array",
            "tickvals": [position for position, _ in visible_ticks],
            "ticktext": [label for _, label in visible_ticks],
        }
    figure = go.Figure(
        go.Heatmap(
            z=transformed.to_numpy(),
            x=[str(value) for value in transformed.columns],
            y=transformed.index.tolist(),
            colorscale=colorscale,
            zmid=zmid,
            zmin=zmin,
            customdata=matrix.to_numpy(),
            colorbar=colorbar,
            hovertemplate=hovertemplate,
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


def section_title_with_info(title: str, description: str) -> None:
    escaped_description = html.escape(description, quote=True)
    st.markdown(
        '<div class="section-title-with-info">'
        f"<h2>{html.escape(title)}</h2>"
        f'<span class="section-info-icon" tabindex="0" aria-label="{escaped_description}">i'
        f'<span class="section-info-tooltip" role="tooltip">{escaped_description}</span>'
        "</span>"
        "</div>",
        unsafe_allow_html=True,
    )


def page_heading(title: str, description: str) -> None:
    st.markdown(
        '<div class="page-heading">'
        f"<h1>{html.escape(title)}</h1>"
        f"<p>{html.escape(description)}</p>"
        "</div>",
        unsafe_allow_html=True,
    )


def render_home() -> None:
    st.markdown(
        """
        # Aedes RNA Atlas

        This site brings together *Aedes aegypti* RNA-seq expression data
        across tissues, feeding conditions, and reproductive states.

        Use the menu above to:

        - **Genes** — search gene symbols and historical identifiers, then compare expression across studies.
        - **Families** — examine IR, OR, GR, and OBP family members.
        - **Compare conditions** — screen genome-wide expression differences between two groups.
        - **Clusters** — inspect relationships between biological samples using PCA, UMAP, or t-SNE.

        ## Data sources

        - [Venkataraman et al., eLife 2023](https://doi.org/10.7554/eLife.80489) — ovary expression across reproductive and drought-resilience states.
        - [Matthews et al., BMC Genomics 2016](https://doi.org/10.1186/s12864-015-2239-0) — female and male tissues across feeding and reproductive conditions.
        - **Nadav Shai · Vosshall lab midgut RNA-seq** — female midgut from non-blood-fed through 72 hours post-blood-meal, plus non-blood-fed male midgut.

        TPM is descriptive normalized abundance. Condition-comparison statistics are
        exploratory; publication-grade differential expression should use raw counts
        with a count-aware model.
        """,
    )


def navigate_home() -> None:
    st.session_state["site_navigation"] = "Home"


navigation_items = ["Home", "Genes", "Families", "Compare conditions", "Clusters"]
if "site_navigation" not in st.session_state:
    st.session_state["site_navigation"] = "Home"
with st.container(key="site_nav"):
    brand_column, navigation_column = st.columns([1.35, 4.65], vertical_alignment="center")
    with brand_column:
        st.button(
            "🧬 Aedes RNA Atlas",
            key="home_logo",
            on_click=navigate_home,
        )
    with navigation_column:
        mode = st.segmented_control(
            "Site navigation",
            navigation_items,
            key="site_navigation",
            label_visibility="collapsed",
            width="stretch",
        )

if mode == "Home":
    render_home()

elif mode == "Genes":
    page_heading(
        "Gene explorer",
        "Search genes and historical identifiers, then compare all resolved genes across the selected experiments.",
    )
    default_selected_keys = [
        key for key in ("elife", "neuro_ru") if key in study_keys
    ] or study_keys[:1]
    selected_keys_for_tokens = [
        key
        for key in st.session_state.get("gene_studies", default_selected_keys)
        if key in study_keys
    ]
    with st.container(key="gene_setup_panel"):
        queries = render_gene_query_editor(
            state_key="gene_query_text",
            editor_key="gene_query_editor",
            selected_keys=selected_keys_for_tokens,
            default_value="ir25a orco",
        )
        selected_keys = st.multiselect(
            "Studies",
            options=study_keys,
            default=default_selected_keys,
            key="gene_studies",
            format_func=lambda key: datasets[key].label,
            help="Choose one or more RNA-seq datasets.",
        )
        scale_column, sort_column, median_column, guide_column = st.columns(4)
        with scale_column:
            x_axis_scale = st.selectbox(
                "TPM scale",
                GENE_X_SCALES,
                index=1,
                help="Controls both the point-plot x-axes and comparison-heatmap color scales. Log scales keep zero values visible using TPM + 1, while ticks and tooltips show the original TPM values.",
            )
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
                help="Show diamonds at each gene's group median.",
            )
        with guide_column:
            show_guides = st.toggle(
                "Show row guides",
                value=True,
                help="Show a dotted horizontal guide on every condition row.",
            )
    if not queries:
        st.info("Enter a gene, for example `Ir25a`, `Orco`, or `AAEL005776`.")
    elif not selected_keys:
        st.info("Choose at least one study.")
    else:
        resolved_for_comparison: dict[str, pd.DataFrame] = {}
        resolved_queries: list[tuple[str, dict[str, pd.DataFrame]]] = []
        combined_summaries: list[pd.DataFrame] = []
        combined_raw_rows: list[pd.DataFrame] = []

        for query in queries:
            per_study: dict[str, pd.DataFrame] = {}
            for key in selected_keys:
                matches = resolve_one(datasets[key], query)
                if not matches.empty:
                    per_study[key] = matches.drop_duplicates("row_id")
                    prior = resolved_for_comparison.get(key)
                    resolved_for_comparison[key] = (
                        per_study[key]
                        if prior is None
                        else pd.concat([prior, per_study[key]], ignore_index=True).drop_duplicates("row_id")
                    )
            resolved_queries.append((query, per_study))

        if resolved_for_comparison:
            render_matched_gene_list(
                matched_gene_entries(resolved_for_comparison),
                len(selected_keys),
            )

        for query, per_study in resolved_queries:
            if not per_study:
                for key in selected_keys:
                    st.warning(
                        f"Gene not found: `{query}` is not present in {datasets[key].label}."
                    )
                continue
            for key in selected_keys:
                if key not in per_study:
                    st.warning(
                        f"Gene not found: `{query}` is not present in {datasets[key].label}."
                    )

        if resolved_for_comparison:
            section_title_with_info(
                "Expression across selected genes",
                "One replicate plot per study. Colors identify genes; diamonds show each gene's group median.",
            )
            for key in selected_keys:
                combined_genes = resolved_for_comparison.get(key)
                if combined_genes is None or combined_genes.empty:
                    continue
                dataset = datasets[key]
                field, field_label = default_grouping(dataset)
                long = expression_long(dataset, combined_genes)
                st.markdown(f"**{dataset.label}**")
                st.plotly_chart(
                    replicate_figure(
                        long,
                        field,
                        field_label,
                        sort_conditions,
                        show_medians,
                        show_guides,
                        x_axis_scale,
                    ),
                    width="stretch",
                    key=f"gene_plot_combined_{key}",
                )

        for query, per_study in resolved_queries:
            if not per_study:
                continue

            for key, matches in per_study.items():
                dataset = datasets[key]
                summary = gene_statistics(dataset, matches)
                summary.insert(0, "Study", dataset.label)
                summary["Alternative names"] = [
                    ", ".join(
                        alias
                        for alias in gene.aliases
                        if alias.casefold() != gene.display_name.casefold()
                    )
                    or "none"
                    for gene in matches.itertuples(index=False)
                ]
                summary["Samples ≥1 TPM"] = summary.apply(
                    lambda row: f"{round(row['detected_pct'] * len(dataset.sample_columns) / 100):.0f}/{len(dataset.sample_columns)}",
                    axis=1,
                )
                combined_summaries.append(summary)
                long = expression_long(dataset, matches)
                long.insert(0, "Study", dataset.label)
                combined_raw_rows.append(long)

        if combined_summaries:
            st.markdown("## Selected gene details")
            summary_table = pd.concat(combined_summaries, ignore_index=True).rename(
                columns={
                    "gene": "Gene",
                    "stable_id": "Stable ID",
                    "raw_symbol": "Published symbol",
                    "family": "Family",
                    "mean_tpm": "Mean TPM",
                    "median_tpm": "Median TPM",
                    "max_tpm": "Maximum TPM",
                    "top_context": "Top context",
                }
            ).drop_duplicates(["Study", "Gene", "Stable ID"])
            st.dataframe(
                summary_table[
                    [
                        "Study",
                        "Gene",
                        "Alternative names",
                        "Family",
                        "Stable ID",
                        "Published symbol",
                        "Median TPM",
                        "Maximum TPM",
                        "Samples ≥1 TPM",
                        "Top context",
                    ]
                ],
                hide_index=True,
                width="stretch",
                column_config={
                    "Median TPM": st.column_config.NumberColumn(format="%.2f"),
                    "Maximum TPM": st.column_config.NumberColumn(format="%.2f"),
                },
            )

            with st.expander("Raw values & download"):
                raw = pd.concat(combined_raw_rows, ignore_index=True)
                raw_table = raw[
                    ["Study", "gene", "stable_id", "sample", "tpm", "tissue", "condition_label", "reproductive_state"]
                ].drop_duplicates(
                    ["Study", "gene", "stable_id", "sample"]
                ).sort_values(["Study", "gene", "tpm"], ascending=[True, True, False]).rename(
                    columns={
                        "gene": "Gene",
                        "stable_id": "Stable ID",
                        "sample": "Sample",
                        "tpm": "TPM",
                        "tissue": "Tissue",
                        "condition_label": "Condition",
                        "reproductive_state": "Reproductive state",
                    }
                )
                st.dataframe(raw_table, hide_index=True, width="stretch", height=280)
                download_tsv(
                    "Download all selected values",
                    raw_table,
                    "selected_genes_expression.tsv",
                    "gene_download_combined",
                )

        if selected_keys and resolved_for_comparison:
            comparisons = []
            for comparison_key in selected_keys:
                comparison_genes = resolved_for_comparison.get(comparison_key)
                if comparison_genes is None or comparison_genes.empty:
                    continue
                dataset = datasets[comparison_key]
                field, _ = default_grouping(dataset)
                long = expression_long(dataset, comparison_genes)
                grouped = grouped_median(long, field)
                comparisons.append((comparison_key, dataset, field, grouped))

            if comparisons:
                section_title_with_info(
                    "Expression heatmaps",
                    "One heatmap per study. The selected TPM scale also controls the colorbars. Colors are scaled within each panel; compare patterns rather than color intensity across papers.",
                )
                for comparison_key, dataset, field, grouped in comparisons:
                    st.markdown(f"**{dataset.label}**")
                    st.plotly_chart(
                        heatmap_figure(
                            grouped,
                            field,
                            "",
                            row_zscore=False,
                            value_scale=x_axis_scale,
                        ),
                        width="stretch",
                        key=f"gene_comparison_heatmap_{comparison_key}",
                    )

elif mode == "Families":
    page_heading(
        "Gene families",
        "Survey chemosensory receptor families across tissues and biological conditions.",
    )
    default_family_keys = [
        key for key in ("neuro_ru", "elife") if key in datasets
    ] or study_keys[:1]
    with st.container(key="family_setup_panel"):
        family_column, studies_column = st.columns(2)
        with family_column:
            family_label = st.selectbox(
                "Gene family",
                [*FAMILIES, CUSTOM_FAMILY_LABEL],
                index=0,
            )
        with studies_column:
            family_keys = st.multiselect(
                "Studies",
                options=study_keys,
                default=default_family_keys,
                key="family_studies",
                format_func=lambda key: datasets[key].label,
            )

        custom_family = family_label == CUSTOM_FAMILY_LABEL
        family_name = FAMILIES.get(family_label, CUSTOM_FAMILY_LABEL)
        family_queries = (
            render_gene_query_editor(
                state_key="family_query_text",
                editor_key="family_query_editor",
                selected_keys=family_keys,
            )
            if custom_family
            else []
        )
        members_by_study = {
            key: (
                resolve_queries(datasets[key], family_queries)
                if custom_family
                else family_members(datasets[key], family_name)
            )
            for key in family_keys
        }
        maximum_family_genes = max(
            (members["display_name"].nunique() for members in members_by_study.values()),
            default=0,
        )

        threshold = 1.0
        count_column, pattern_column = st.columns(2)
        with count_column:
            if maximum_family_genes:
                count_options = gene_count_options(maximum_family_genes)
                top_n = st.select_slider(
                    "Genes in each heatmap",
                    options=count_options,
                    value=maximum_family_genes,
                    format_func=lambda value: (
                        f"All ({maximum_family_genes})"
                        if value == maximum_family_genes
                        else str(value)
                    ),
                    help="For each study, genes are ranked by mean TPM across all biological samples, the standard summary for overall expression. Conditions with more replicates therefore contribute more samples. Selecting N shows exactly the N highest-ranked genes; All shows every matched gene. Heatmap cells still show the median TPM within each condition.",
                )
            else:
                top_n = 0
                st.select_slider(
                    "Genes in each heatmap",
                    options=[0],
                    value=0,
                    format_func=lambda _: "All (0)",
                    disabled=True,
                    help="Add genes to the custom family to enable this control.",
                )
        with pattern_column:
            row_zscore = st.toggle(
                "Show relative pattern within each gene",
                value=True,
                help="Convert each gene's log₂(TPM + 1) values to z-scores across conditions. Positive means above that gene's average; negative means below it.",
            )

    if custom_family and family_queries:
        matched_custom_genes = {
            key: members
            for key, members in members_by_study.items()
            if not members.empty
        }
        render_matched_gene_list(
            matched_gene_entries(matched_custom_genes),
            len(family_keys),
        )
        for query in family_queries:
            for key in family_keys:
                if resolve_one(datasets[key], query).empty:
                    st.warning(
                        f"Gene not found: `{query}` is not present in {datasets[key].label}."
                    )

    if not family_keys:
        st.info("Choose at least one study.")
    elif custom_family and not family_queries:
        st.info("Enter genes or identifiers to create a custom family.")
    else:
        coverage = []
        family_data: dict[str, tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, str]] = {}
        for key in family_keys:
            dataset = datasets[key]
            members = members_by_study[key]
            if members.empty:
                coverage.append({"Study": dataset.label, "Family genes": 0, f"Samples ≥{threshold:g} TPM": "—"})
                continue
            family_gene_count = members["display_name"].nunique()
            field, _ = default_grouping(dataset)
            long = expression_long(dataset, members)
            group_medians = grouped_median(long, field)
            mean_expression = (
                long.groupby("gene", as_index=False)["tpm"].mean().rename(
                    columns={"tpm": "mean_tpm"}
                )
            )
            detected = (
                long.groupby(["gene", "sample"], as_index=False)["tpm"]
                .max()
                .assign(detected=lambda frame: frame["tpm"] >= threshold)
                .groupby("gene", as_index=False)["detected"]
                .sum()
                .rename(columns={"detected": "Detected n"})
            )
            ranking = (
                gene_statistics(dataset, members)
                .sort_values("mean_tpm", ascending=False)
                .drop_duplicates("gene")
                .drop(columns="mean_tpm")
                .merge(mean_expression, on="gene", how="left")
                .merge(detected, on="gene", how="left")
            )
            detection_column = f"Samples ≥{threshold:g} TPM"
            ranking[detection_column] = (
                ranking["Detected n"].astype(int).astype(str)
                + f"/{len(dataset.sample_columns)}"
            )
            ranking = ranking.sort_values(
                ["mean_tpm", "gene"],
                ascending=[False, True],
            )

            selected_names = ranking["gene"].head(top_n).tolist()
            selected_members = members[members["display_name"].isin(selected_names)]
            selected_long = expression_long(dataset, selected_members)
            selected_grouped = grouped_median(selected_long, field)
            selected_grouped["gene"] = pd.Categorical(
                selected_grouped["gene"],
                categories=selected_names,
                ordered=True,
            )
            family_data[key] = (members, ranking, selected_grouped, field)
            coverage.append(
                {
                    "Study": dataset.label,
                    "Family genes": family_gene_count,
                    f"Samples ≥{threshold:g} TPM": f"{(ranking['Detected n'] > 0).sum()}/{family_gene_count} genes",
                }
            )

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
                    "mean_tpm": "Mean TPM",
                    "top_context": "Top context",
                }
            )[
                ["Gene", "Stable ID", "Mean TPM", detection_column, "Top context"]
            ]
            st.dataframe(
                concise.head(top_n),
                hide_index=True,
                width="stretch",
                column_config={
                    "Mean TPM": st.column_config.NumberColumn(format="%.2f")
                },
            )
            with st.expander("Complete family table & download"):
                st.dataframe(concise, hide_index=True, width="stretch", height=340)
                download_tsv(
                    "Download complete family matrix",
                    matrix_for_genes(dataset, members),
                    f"{key}_{'custom' if custom_family else family_name.split('(')[-1].rstrip(')').casefold()}_family_tpm.tsv",
                    f"family_download_{key}",
                )

elif mode == "Compare conditions":
    page_heading(
        "Compare conditions",
        "Screen every measured gene for expression differences between two biological groups.",
    )
    st.caption(
        "This is an MA plot: right means higher average TPM. A ratio of 1× means equal expression; above 1× means higher in A; below 1× means higher in B."
    )
    with st.container(key="condition_setup_panel"):
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
        filter_text = st.text_input(
            "Filter results by gene or Stable ID",
            placeholder="e.g. Ir25a or AAEL005776",
        ).strip()

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

else:
    page_heading(
        "Sample clusters",
        "Map whole-transcriptome similarity between biological samples with PCA, UMAP, or t-SNE.",
    )
    st.caption(
        "Each point is one biological sample. Nearby points have more similar whole-transcriptome expression profiles."
    )
    with st.container(key="cluster_setup_panel"):
        study_column, method_column = st.columns(2)
        with study_column:
            cluster_key = st.selectbox(
                "Study",
                options=study_keys,
                format_func=lambda key: datasets[key].label,
                key="cluster_study",
            )
        with method_column:
            cluster_method = st.selectbox("Method", METHODS, key="cluster_method")

        cluster_dataset = datasets[cluster_key]
        if cluster_key.startswith("neuro_"):
            color_candidates = [
                ("tissue", "Tissue"),
                ("condition_label", "Feeding / reproductive state"),
                ("sex", "Sex"),
                ("tissue_condition", "Tissue + condition"),
            ]
        else:
            color_candidates = [
                ("reproductive_state", "Reproductive state"),
            ]
        color_options = [
            field
            for field, _ in color_candidates
            if field in cluster_dataset.samples
            and cluster_dataset.samples[field].fillna("").astype(str).nunique() > 1
        ]
        color_labels = dict(color_candidates)
        if not color_options:
            color_options = ["sample"]
            color_labels["sample"] = "Sample"

        color_column, genes_column = st.columns(2)
        with color_column:
            color_field = st.selectbox(
                "Color by",
                options=color_options,
                format_func=lambda field: color_labels[field],
                key="cluster_color",
            )
        with genes_column:
            variable_genes = st.slider(
                "Most-variable genes",
                min_value=250,
                max_value=5_000,
                value=2_000,
                step=250,
                help="Genes are ranked by variance after log-transforming TPM. Constant genes are excluded.",
            )

    try:
        cluster_frame, cluster_x, cluster_y, cluster_details = cluster_embedding_cached(
            DATA_SCHEMA_VERSION,
            cluster_key,
            cluster_method,
            variable_genes,
        )
    except ValueError as exc:
        st.warning(str(exc))
    else:
        cluster_frame[color_field] = (
            cluster_frame[color_field].fillna("").astype(str).replace("", "Unspecified")
        )
        color_order = cluster_frame[color_field].drop_duplicates().tolist()
        hover_fields = [
            field
            for field in ("sample", "tissue", "condition_label", "sex")
            if field in cluster_frame.columns
        ]
        hover_data = {field: True for field in hover_fields}
        hover_data.update({"x": False, "y": False})
        cluster_figure = px.scatter(
            cluster_frame,
            x="x",
            y="y",
            color=color_field,
            hover_name="sample",
            hover_data=hover_data,
            labels={
                "x": cluster_x,
                "y": cluster_y,
                color_field: color_labels[color_field],
            },
            category_orders={color_field: color_order},
            color_discrete_sequence=px.colors.qualitative.Set2
            + px.colors.qualitative.Pastel,
        )
        cluster_figure.update_traces(
            marker={"size": 11, "opacity": 0.9, "line": {"width": 0.6, "color": "#111827"}}
        )
        cluster_figure.update_layout(
            height=650,
            margin={"l": 30, "r": 25, "t": 25, "b": 55},
            legend={"title": {"text": color_labels[color_field]}},
        )
        cluster_figure.update_yaxes(showgrid=False, zeroline=False)
        section_title_with_info(
            "Cluster map",
            f"{cluster_method} · {cluster_details}. TPM was transformed as log₂(TPM + 1), then each selected gene was standardized across samples. PCA preserves broad linear variation. UMAP and t-SNE emphasize local neighborhoods; their axis values and distances between far-apart groups are not directly interpretable.",
        )
        st.plotly_chart(
            cluster_figure,
            width="stretch",
            key=f"cluster_plot_{cluster_key}_{cluster_method}_{variable_genes}_{color_field}",
        )
