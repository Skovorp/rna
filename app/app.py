from __future__ import annotations

import hashlib
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
from expression_explorer.data import (
    DATASET_ORDER,
    expression_long,
    family_members,
    gene_statistics,
    load_datasets,
    matrix_for_genes,
    search_genes,
)
from expression_explorer.differential import (
    condition_label,
    contrast_label,
    contrast_sample_counts,
    load_differential_contrasts,
    load_differential_results,
)


APP_DIR = Path(__file__).resolve().parent
EXPRESSION_DIR = APP_DIR.parent / "expression"
DATA_SCHEMA_VERSION = "2026-08-04-midgut-deseq2-v1"

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

SAMPLE_GROUP_COLOR_SEQUENCE = px.colors.qualitative.Set2 + px.colors.qualitative.Pastel

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
    .st-key-differential_setup_panel,
    .st-key-cluster_setup_panel {
        margin: .7rem 0 1rem;
        padding: .9rem 1rem .7rem;
        border: 1px solid rgba(148, 163, 184, .25);
        border-radius: .75rem;
        background: rgba(148, 163, 184, .09);
    }
    .matched-genes-title {
        margin-top: .35rem;
        font-size: .82rem;
        font-weight: 750;
    }
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
differential_contrasts = load_differential_contrasts(EXPRESSION_DIR)
differential_study_keys = sorted(
    study_keys,
    key=lambda key: (key not in differential_contrasts, study_keys.index(key)),
)


@st.cache_data(show_spinner="Computing sample map…")
def cluster_embedding_cached(
    schema_version: str,
    dataset_key: str,
    method: str,
    variable_genes: int,
):
    dataset = datasets_resource(schema_version)[dataset_key]
    return sample_embedding(dataset, method, variable_genes)


@st.cache_data(show_spinner="Loading differential-expression results…")
def differential_results_cached(
    schema_version: str,
    dataset_key: str,
    contrast_id: str,
):
    dataset = datasets_resource(schema_version)[dataset_key]
    contrast = next(
        contrast
        for contrast in differential_contrasts[dataset_key]
        if contrast.contrast_id == contrast_id
    )
    return load_differential_results(dataset, contrast)


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


def mean_expression_by_study(
    per_study: dict[str, pd.DataFrame],
    selected_study_keys: list[str],
) -> dict[str, dict[str, float]]:
    """Return one mean-TPM mapping per selected study."""
    means_by_study: dict[str, dict[str, float]] = {}
    for study_key in selected_study_keys:
        matches = per_study.get(study_key)
        if matches is None or matches.empty:
            means_by_study[study_key] = {}
            continue
        long = expression_long(datasets[study_key], matches)
        means_by_study[study_key] = {
            str(gene_name).casefold(): float(mean_tpm)
            for gene_name, mean_tpm in long.groupby("gene")["tpm"].mean().items()
        }
    return means_by_study


def render_matched_gene_table(
    entries: list[dict[str, object]],
    selected_study_count: int,
    state_key: str,
    mean_tpm_by_study: dict[str, dict[str, float]],
) -> set[str]:
    if not entries:
        return set()

    enabled_by_gene = st.session_state.setdefault(state_key, {})
    entry_keys = [str(entry["display_name"]).casefold() for entry in entries]
    for entry_key in entry_keys:
        enabled_by_gene.setdefault(entry_key, True)

    revision_key = f"{state_key}_revision"
    st.session_state.setdefault(revision_key, 0)
    all_on_column, all_off_column, _ = st.columns([1, 1, 5])
    with all_on_column:
        if st.button("Turn all on", key=f"{state_key}_all_on", width="stretch"):
            enabled_by_gene.update(dict.fromkeys(entry_keys, True))
            st.session_state[revision_key] += 1
    with all_off_column:
        if st.button("Turn all off", key=f"{state_key}_all_off", width="stretch"):
            enabled_by_gene.update(dict.fromkeys(entry_keys, False))
            st.session_state[revision_key] += 1

    mean_columns = {
        study_key: f"Mean TPM · {datasets[study_key].label}"
        for study_key in mean_tpm_by_study
    }
    rows: list[dict[str, object]] = []
    for entry in entries:
        display_name = str(entry["display_name"])
        display_key = display_name.casefold()
        families = list(entry["families"].values())
        aliases = sorted(entry["aliases"].values(), key=str.casefold)
        found_count = len(entry["study_keys"])
        row: dict[str, object] = {
            "Include": bool(enabled_by_gene[display_key]),
            "Gene": display_name,
            "Family": "/".join(families),
            "Alternative names": (
                ", ".join(aliases) if aliases else "No alternative names found"
            ),
            "Study coverage": f"{found_count}/{selected_study_count} studies",
        }
        for study_key, column_name in mean_columns.items():
            row[column_name] = mean_tpm_by_study[study_key].get(display_key, np.nan)
        rows.append(row)

    first_mean_column = next(iter(mean_columns.values()))
    rows.sort(
        key=lambda row: (
            pd.isna(row[first_mean_column]),
            -float(row[first_mean_column])
            if not pd.isna(row[first_mean_column])
            else 0.0,
            str(row["Gene"]).casefold(),
        )
    )
    entry_keys = [str(row["Gene"]).casefold() for row in rows]
    st.markdown(
        f'<div class="matched-genes-title">Matched genes · {len(entries)}</div>',
        unsafe_allow_html=True,
    )
    table = pd.DataFrame(rows)
    row_fingerprint = hashlib.sha1(
        "\n".join(entry_keys).encode("utf-8")
    ).hexdigest()[:12]
    edited = st.data_editor(
        table,
        hide_index=True,
        width="stretch",
        height=min(38 + 35 * len(table), 318),
        row_height=35,
        disabled=[
            "Gene",
            "Family",
            *mean_columns.values(),
            "Alternative names",
            "Study coverage",
        ],
        column_config={
            "Include": st.column_config.CheckboxColumn(
                "Include",
                help="Turn a gene off to remove it from every result below.",
                width="small",
            ),
            "Gene": st.column_config.TextColumn(width="small"),
            "Family": st.column_config.TextColumn(width="small"),
            **{
                column_name: st.column_config.NumberColumn(
                    format="%.2f",
                    help=f"Mean TPM across every biological sample in {datasets[study_key].label}.",
                    width="small",
                )
                for study_key, column_name in mean_columns.items()
            },
            "Alternative names": st.column_config.TextColumn(width="large"),
            "Study coverage": st.column_config.TextColumn(width="small"),
        },
        key=(
            f"{state_key}_editor_{row_fingerprint}_"
            f"{st.session_state[revision_key]}"
        ),
    )
    enabled_by_gene.update(
        {
            str(row.Gene).casefold(): bool(row.Include)
            for row in edited.itertuples(index=False)
        }
    )
    enabled_count = sum(enabled_by_gene[entry_key] for entry_key in entry_keys)
    st.caption(
        f"{enabled_count} of {len(entries)} matched genes included in the results below."
    )
    return {
        entry_key
        for entry_key in entry_keys
        if enabled_by_gene[entry_key]
    }


def filter_matched_genes(
    per_study: dict[str, pd.DataFrame],
    enabled_gene_keys: set[str],
) -> dict[str, pd.DataFrame]:
    filtered: dict[str, pd.DataFrame] = {}
    for study_key, matches in per_study.items():
        included = matches[
            matches["display_name"].astype(str).str.casefold().isin(enabled_gene_keys)
        ]
        if not included.empty:
            filtered[study_key] = included
    return filtered


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


def cluster_gene_count_options(maximum: int) -> list[int]:
    """Return useful clustering gene-count stops followed by the complete set."""
    standard_options = (250, 500, 1_000, 2_000, 5_000, 10_000)
    return [value for value in standard_options if value < maximum] + [maximum]


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


def sample_groupings(dataset) -> list[tuple[str, str]]:
    """Return the biological sample metadata exposed for filtering and color."""
    configured = {
        "elife": [
            ("reproductive_state", "Reproductive state / time"),
        ],
        "neuro_ru": [
            ("sex", "Sex"),
            ("tissue", "Anatomy"),
            ("condition_label", "Feeding / reproductive state"),
        ],
        "midgut": [
            ("sex", "Sex"),
            ("timepoint", "Blood-meal time"),
            ("condition_label", "Sex + blood-meal time"),
        ],
    }.get(dataset.key, [])
    return [
        (field, label)
        for field, label in configured
        if field in dataset.samples.columns
        and dataset.samples[field].fillna("").astype(str).nunique() > 1
    ]


def sample_group_values(dataset, field: str) -> list[str]:
    values = dataset.samples[field].fillna("").astype(str).replace("", "Unspecified")
    return values.drop_duplicates().tolist()


def render_gene_sample_controls(dataset) -> tuple[dict[str, list[str]], str | None, str | None]:
    groupings = sample_groupings(dataset)
    if not groupings:
        return {}, None, None

    filters: dict[str, list[str]] = {}
    labels = dict(groupings)
    columns = st.columns(len(groupings) + 1)
    for column, (field, label) in zip(columns, groupings):
        with column:
            filters[field] = st.multiselect(
                f"Filter · {label}",
                options=sample_group_values(dataset, field),
                default=[],
                key=f"gene_filter_{dataset.key}_{field}",
                help="Leave empty to include every value in this group.",
            )
    with columns[-1]:
        color_field = st.selectbox(
            "Color condition labels by",
            options=[None, *labels],
            index=0,
            key=f"gene_label_color_{dataset.key}",
            format_func=lambda field: "None" if field is None else labels[field],
            help="None leaves every condition label in the default text color. Point colors always identify genes.",
        )
    return filters, color_field, labels.get(color_field)


def filter_expression_samples(
    long: pd.DataFrame, filters: dict[str, list[str]]
) -> pd.DataFrame:
    filtered = long.copy()
    for field, selected_values in filters.items():
        if not selected_values:
            continue
        values = filtered[field].fillna("").astype(str).replace("", "Unspecified")
        filtered = filtered.loc[values.isin(selected_values)].copy()
    return filtered


def category_label_colors(
    long: pd.DataFrame,
    category_field: str,
    color_field: str | None,
) -> dict[str, str]:
    """Map plot categories to colors from one biological metadata field."""
    if not color_field:
        return {}
    color_values = long[color_field].fillna("").astype(str).replace("", "Unspecified")
    color_order = color_values.drop_duplicates().tolist()
    color_map = {
        value: SAMPLE_GROUP_COLOR_SEQUENCE[index % len(SAMPLE_GROUP_COLOR_SEQUENCE)]
        for index, value in enumerate(color_order)
    }
    categories = long[category_field].fillna("Unspecified").astype(str)
    result: dict[str, str] = {}
    for category in categories.drop_duplicates():
        groups = color_values.loc[categories.eq(category)].drop_duplicates().tolist()
        result[category] = color_map[groups[0]] if len(groups) == 1 else "#9aa8a5"
    return result


def category_order_by_group(
    long: pd.DataFrame,
    category_field: str,
    color_field: str | None,
) -> list[str]:
    """Group categories by selected metadata, preserving original order within groups."""
    categories = long[category_field].fillna("Unspecified").astype(str)
    original_order = categories.drop_duplicates().tolist()
    if not color_field:
        return original_order

    color_values = long[color_field].fillna("").astype(str).replace("", "Unspecified")
    group_order = color_values.drop_duplicates().tolist()
    group_rank = {value: index for index, value in enumerate(group_order)}
    original_rank = {value: index for index, value in enumerate(original_order)}

    def category_group_rank(category: str) -> int:
        groups = color_values.loc[categories.eq(category)].drop_duplicates().tolist()
        return min((group_rank[group] for group in groups), default=len(group_rank))

    return sorted(
        original_order,
        key=lambda category: (category_group_rank(category), original_rank[category]),
    )


def grouped_median(
    long: pd.DataFrame,
    field: str,
    order: list[str] | None = None,
) -> pd.DataFrame:
    if order is None:
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
    show_medians: bool = True,
    show_guides: bool = True,
    x_scale: str = "Log base 2",
    category_colors: dict[str, str] | None = None,
    condition_order: list[str] | None = None,
) -> go.Figure:
    plot = long.copy()
    plot["axis_tpm"] = tpm_axis_position(plot["tpm"], x_scale)
    plot[field] = plot[field].fillna("Unspecified").astype(str)
    if condition_order is None:
        condition_order = plot[field].drop_duplicates().tolist()
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
        legend={
            "title": {"text": ""},
            "itemclick": False,
            "itemdoubleclick": False,
        },
    )
    if category_colors:
        figure.update_yaxes(
            tickmode="array",
            tickvals=condition_order,
            ticktext=condition_order,
            tickfont={"color": "rgba(0,0,0,0)"},
        )
        for condition in condition_order:
            figure.add_annotation(
                x=0,
                xref="paper",
                xshift=-8,
                y=condition,
                yref="y",
                text=html.escape(str(condition)),
                showarrow=False,
                xanchor="right",
                yanchor="middle",
                align="right",
                font={"color": category_colors.get(condition, "#9aa8a5"), "size": 12},
            )
    return figure


def heatmap_figure(
    grouped: pd.DataFrame,
    field: str,
    title: str,
    row_zscore: bool,
    value_scale: str = "Log base 2",
    category_colors: dict[str, str] | None = None,
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
    bottom_margin = 125 if category_colors else 60
    figure.update_layout(
        title=title,
        height=max(390, min(1050, 145 + 23 * len(transformed))),
        margin={"l": 20, "r": 25, "t": 55, "b": bottom_margin},
        xaxis={"tickangle": -35},
        yaxis={"autorange": "reversed"},
    )
    if category_colors:
        categories = [str(value) for value in transformed.columns]
        figure.update_xaxes(
            tickmode="array",
            tickvals=categories,
            ticktext=categories,
            tickfont={"color": "rgba(0,0,0,0)"},
        )
        for category in categories:
            figure.add_annotation(
                x=category,
                xref="x",
                y=0,
                yref="paper",
                yshift=-10,
                text=html.escape(category),
                textangle=-35,
                showarrow=False,
                xanchor="right",
                yanchor="top",
                align="right",
                font={"color": category_colors.get(category, "#9aa8a5"), "size": 12},
            )
    return figure


def differential_log2_fold_change_range(results: pd.DataFrame) -> list[float]:
    fold_changes = results.loc[
        results["ma_plot_eligible"], "log2_fold_change"
    ].dropna()
    if fold_changes.empty:
        return [-2.0, 2.0]
    limit = max(2.0, float(np.ceil(fold_changes.abs().quantile(0.995))))
    return [-limit, limit]


def differential_abundance_range(results: pd.DataFrame) -> list[float]:
    finite_abundance = results.loc[
        results["ma_plot_eligible"], "base_mean"
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


def differential_ma_figure(results: pd.DataFrame, fdr_threshold: float) -> go.Figure:
    plotted = results[results["ma_plot_eligible"]].copy()
    fold_change_range = differential_log2_fold_change_range(results)
    abundance_range = differential_abundance_range(results)
    abundance_ticks, abundance_labels = base10_ticks(abundance_range)
    plotted["passes_fdr"] = plotted["fdr"] < fdr_threshold
    figure = go.Figure()
    for passes_fdr, label, color in (
        (False, f"FDR ≥ {fdr_threshold:g}", "#66706f"),
        (True, f"FDR < {fdr_threshold:g}", "#f5b85b"),
    ):
        subset = plotted[plotted["passes_fdr"].eq(passes_fdr)]
        figure.add_trace(
            go.Scattergl(
                x=subset["base_mean"],
                y=subset["log2_fold_change"],
                mode="markers",
                name=label,
                marker={"size": 3, "color": color, "opacity": 1.0},
                customdata=subset[
                    [
                        "gene",
                        "stable_id",
                        "base_mean",
                        "fold_change",
                        "lfc_se",
                        "fdr",
                    ]
                ].to_numpy(),
                hovertemplate=(
                    "<b>%{customdata[0]}</b><br>"
                    "Stable ID: %{customdata[1]}<br>"
                    "DESeq2 base mean: %{customdata[2]:.3f}<br>"
                    "Log₂ fold change: %{y:.3f}<br>"
                    "Fold change (target / reference): %{customdata[3]:.3f}×<br>"
                    "Log₂ fold-change SE: %{customdata[4]:.3f}<br>"
                    "FDR: %{customdata[5]:.3g}<extra></extra>"
                ),
            )
        )
    figure.add_hline(
        y=0,
        line={"color": "rgba(148,163,184,.38)", "dash": "dot", "width": 1},
    )
    figure.update_layout(
        height=510,
        margin={"l": 25, "r": 25, "t": 35, "b": 55},
        xaxis={
            "title": "DESeq2 base mean (logarithmic scale)",
            "type": "log",
            "range": abundance_range,
            "tickmode": "array",
            "tickvals": abundance_ticks,
            "ticktext": abundance_labels,
        },
        yaxis={
            "title": "Log₂ fold change · target / reference",
            "range": fold_change_range,
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
        - **Differential expression** — browse precomputed DESeq2 contrasts from nf-core/differentialabundance.
        - **Clusters** — inspect relationships between biological samples using PCA, UMAP, or t-SNE.

        ## Data sources

        - [Venkataraman et al., eLife 2023](https://doi.org/10.7554/eLife.80489) — ovary expression across reproductive and drought-resilience states.
        - [Matthews et al., BMC Genomics 2016](https://doi.org/10.1186/s12864-015-2239-0) — female and male tissues across feeding and reproductive conditions.
        - **Nadav Shai · Vosshall lab midgut RNA-seq** — female midgut from non-blood-fed through 72 hours post-blood-meal, plus non-blood-fed male midgut.

        ## Reference, annotation, and provenance

        Our raw-FASTQ analyses use one shared *Aedes aegypti* reference pair:

        - **Genome FASTA:** `VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta`
          (AaegL5 / AaegLVP_AGWG; SHA-256
          `fd96bcf4d05fa54b0dc4edeebaf077ab6206c35bda1faa150a392dfdcd0545ec`).
        - **Gene annotation GTF:**
          `AaegLVP_VB58-Jove19_MT_noS1_geneNames.sorted.gtf`
          (VectorBase 58 with Jové et al. 2019 gene names; SHA-256
          `0bf20c3fae7f8788e44b56fdbc1e81f5dd502da8c4350d5e4c8757a048a74580`).

        The ovary FASTQs were quantified by us with **Salmon 2.4.1**. The
        midgut run used Trim Galore followed by **Salmon 1.10.3
        pseudoalignment and quantification**. It did not create genomic BAM
        alignments: the recorded nf-core setting was `skip_alignment: true`
        with `pseudo_aligner: salmon`.

        | Dataset | What this site currently displays | Processing provenance |
        |---|---|---|
        | Venkataraman et al. 2023 · ovary (`PRJNA796320`) | Our Salmon 2.4.1 gene-level TPM matrix | We quantified all 33 biological samples from the raw paired-end FASTQs against the shared AaegL5 FASTA and VectorBase 58 + Jové GTF above. The study's published TPM values are not displayed. |
        | Matthews et al. 2016 · neurotranscriptome (`PRJNA236239`) | Published `AaegL.RU` TPM matrix from the study supplement | The displayed values retain the authors' annotation. Our raw-FASTQ reprocessing is tracked separately and is not presented as complete here. |
        | Nadav Shai · midgut | Our Salmon gene-level TPM matrix and our seven DESeq2 contrasts | We ran nf-core/rnaseq 3.26.0 and nf-core/differentialabundance 2.0.0 ourselves from the raw paired-end FASTQs using the reference pair above. |

        TPM is descriptive normalized abundance. Differential-expression statistics
        are displayed only when precomputed count-aware nf-core results are available.

        ## Exact nf-core commands

        These are the original commands and worker paths preserved from the
        completed midgut run. Nextflow was pinned to 25.10.7. The referenced
        parameter files are reproduced directly below each command.
        """,
    )

    st.markdown("### RNA-seq quantification · nf-core/rnaseq 3.26.0")
    st.code(
        """export PATH="/opt/conda/bin:/usr/local/bin:$PATH"
export NXF_HOME=/rna/.nextflow
export NXF_VER=25.10.7
export NXF_CONDA_CACHEDIR=/rna/conda
export NXF_TEMP=/rna/tmp
export NXF_OPTS='-Xms1g -Xmx8g'

PROJECT_ROOT=/rna/project/nadav_shai
WORK_ROOT=/rna/work/nadav_shai_nfcore

nextflow run nf-core/rnaseq \\
    -r 3.26.0 \\
    -profile mamba \\
    -c "$PROJECT_ROOT/vast.config" \\
    -params-file "$PROJECT_ROOT/nfcore-params.json" \\
    -work-dir "$WORK_ROOT" \\
    -resume""",
        language="bash",
    )
    with st.expander("Exact nf-core/rnaseq parameter file"):
        st.code(
            """{
  "input": "/rna/project/nadav_shai/samplesheet.csv",
  "fasta": "/rna/input/Nadav_Shai_Midgut_RNAseq/VectorBase-63_AaegyptiLVP_AGWG_Genome/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta",
  "gtf": "/rna/reference/AaegLVP_VB58-Jove19_MT_noS1_geneNames.sorted.gtf",
  "pseudo_aligner": "salmon",
  "skip_alignment": true,
  "save_reference": true,
  "max_cpus": 28,
  "max_memory": "170.GB",
  "max_time": "72.h",
  "outdir": "/rna/results/nadav_shai_nfcore"
}""",
            language="json",
        )

    st.markdown(
        "### Differential expression · nf-core/differentialabundance 2.0.0"
    )
    st.code(
        """export PATH="/opt/conda/bin:/usr/local/bin:$PATH"
export NXF_HOME=/rna/.nextflow
export NXF_VER=25.10.7
export NXF_CONDA_CACHEDIR=/rna/conda
export NXF_TEMP=/rna/tmp
export NXF_OPTS='-Xms1g -Xmx8g'

PROJECT_ROOT=/rna/project/nadav_shai
WORK_ROOT=/rna/work/nadav_shai_differential

nextflow run nf-core/differentialabundance \\
    -r 2.0.0 \\
    -profile rnaseq,mamba \\
    -c "$PROJECT_ROOT/vast.config" \\
    -params-file "$PROJECT_ROOT/differential-params.json" \\
    -work-dir "$WORK_ROOT" \\
    -resume""",
        language="bash",
    )
    with st.expander("Exact differential-expression parameters and contrasts"):
        st.markdown("**`differential-params.json`**")
        st.code(
            """{
  "study_name": "nadav_shai_midgut_rnaseq",
  "study_type": "rnaseq",
  "input": "/rna/project/nadav_shai/observations.csv",
  "contrasts": "/rna/project/nadav_shai/contrasts.csv",
  "matrix": "/rna/results/nadav_shai_nfcore/salmon/salmon.merged.gene_counts.tsv",
  "feature_length_matrix": "/rna/results/nadav_shai_nfcore/salmon/salmon.merged.gene_lengths.tsv",
  "gtf": "/rna/reference/AaegLVP_VB58-Jove19_MT_noS1_geneNames.sorted.gtf",
  "seed": 20260804,
  "outdir": "/rna/results/nadav_shai_differential"
}""",
            language="json",
        )
        st.markdown("**`contrasts.csv`**")
        st.code(
            """id,variable,reference,target
female_3hBF_vs_NBF,condition,female_NBF,female_3hBF
female_6hBF_vs_NBF,condition,female_NBF,female_6hBF
female_12hBF_vs_NBF,condition,female_NBF,female_12hBF
female_24hBF_vs_NBF,condition,female_NBF,female_24hBF
female_48hBF_vs_NBF,condition,female_NBF,female_48hBF
female_72hBF_vs_NBF,condition,female_NBF,female_72hBF
male_vs_female_NBF,condition,female_NBF,male_midgut""",
            language="text",
        )

    st.caption(
        "Recorded software: Salmon 1.10.3; nf-core/rnaseq 3.26.0; "
        "nf-core/differentialabundance 2.0.0; DESeq2 1.34.0; Nextflow 25.10.7."
    )


def navigate_home() -> None:
    st.session_state["site_navigation"] = "Home"


navigation_items = ["Home", "Genes", "Families", "Differential expression", "Clusters"]
requested_mode = st.query_params.get("page", "Home")
if isinstance(requested_mode, list):
    requested_mode = requested_mode[-1]
if requested_mode not in navigation_items:
    requested_mode = "Home"
if "site_navigation" not in st.session_state:
    st.session_state["site_navigation"] = requested_mode
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
if mode in navigation_items and st.query_params.get("page") != mode:
    st.query_params["page"] = mode

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
        scale_column, median_column, guide_column = st.columns(3)
        with scale_column:
            x_axis_scale = st.selectbox(
                "TPM scale",
                GENE_X_SCALES,
                index=1,
                help="Controls both the point-plot x-axes and comparison-heatmap color scales. Log scales keep zero values visible using TPM + 1, while ticks and tooltips show the original TPM values.",
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
        all_resolved_for_comparison: dict[str, pd.DataFrame] = {}
        resolved_queries: list[tuple[str, dict[str, pd.DataFrame]]] = []
        combined_summaries: list[pd.DataFrame] = []
        combined_raw_rows: list[pd.DataFrame] = []

        for query in queries:
            per_study: dict[str, pd.DataFrame] = {}
            for key in selected_keys:
                matches = resolve_one(datasets[key], query)
                if not matches.empty:
                    per_study[key] = matches.drop_duplicates("row_id")
                    prior = all_resolved_for_comparison.get(key)
                    all_resolved_for_comparison[key] = (
                        per_study[key]
                        if prior is None
                        else pd.concat([prior, per_study[key]], ignore_index=True).drop_duplicates("row_id")
                    )
            resolved_queries.append((query, per_study))

        resolved_for_comparison: dict[str, pd.DataFrame] = {}
        if all_resolved_for_comparison:
            enabled_gene_keys = render_matched_gene_table(
                matched_gene_entries(all_resolved_for_comparison),
                len(selected_keys),
                "gene_matched_gene_selection",
                mean_expression_by_study(all_resolved_for_comparison, selected_keys),
            )
            resolved_for_comparison = filter_matched_genes(
                all_resolved_for_comparison,
                enabled_gene_keys,
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

        if all_resolved_for_comparison and not resolved_for_comparison:
            st.info("Turn on at least one matched gene to show expression results.")

        if resolved_for_comparison:
            section_title_with_info(
                "Expression across selected genes",
                "Each study has sample filters, optional condition-label coloring, a replicate plot, and its heatmap. Point colors identify genes; diamonds show each gene's group median.",
            )
            for key in selected_keys:
                combined_genes = resolved_for_comparison.get(key)
                if combined_genes is None or combined_genes.empty:
                    continue
                dataset = datasets[key]
                field, field_label = default_grouping(dataset)
                long = expression_long(dataset, combined_genes)
                st.markdown(f"### {dataset.label}")
                filters, label_color_field, _ = render_gene_sample_controls(dataset)
                filtered_long = filter_expression_samples(long, filters)
                filtered_sample_count = filtered_long["sample"].nunique()
                if any(filters.values()):
                    st.caption(
                        f"Showing {filtered_sample_count} of {len(dataset.sample_columns)} samples."
                    )
                else:
                    st.caption(f"Showing all {len(dataset.sample_columns)} samples.")
                if filtered_long.empty:
                    st.warning("No samples match these filters.")
                    continue
                label_colors = category_label_colors(
                    filtered_long,
                    field,
                    label_color_field,
                )
                condition_order = category_order_by_group(
                    filtered_long,
                    field,
                    label_color_field,
                )
                st.plotly_chart(
                    replicate_figure(
                        filtered_long,
                        field,
                        field_label,
                        show_medians=show_medians,
                        show_guides=show_guides,
                        x_scale=x_axis_scale,
                        category_colors=label_colors,
                        condition_order=condition_order,
                    ),
                    width="stretch",
                    key=f"gene_plot_combined_{key}",
                )
                grouped = grouped_median(filtered_long, field, condition_order)
                st.plotly_chart(
                    heatmap_figure(
                        grouped,
                        field,
                        "",
                        row_zscore=False,
                        value_scale=x_axis_scale,
                        category_colors=label_colors,
                    ),
                    width="stretch",
                    key=f"gene_comparison_heatmap_{key}",
                )

        for query, per_study in resolved_queries:
            if not per_study:
                continue

            for key, matches in per_study.items():
                matches = matches[
                    matches["display_name"]
                    .astype(str)
                    .str.casefold()
                    .isin(enabled_gene_keys)
                ]
                if matches.empty:
                    continue
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
        all_members_by_study = {
            key: (
                resolve_queries(datasets[key], family_queries)
                if custom_family
                else family_members(datasets[key], family_name)
            )
            for key in family_keys
        }
        matched_family_entries = matched_gene_entries(all_members_by_study)
        enabled_family_gene_keys: set[str] = set()
        if matched_family_entries:
            family_selection_id = (
                "custom"
                if custom_family
                else family_label.split(" · ", maxsplit=1)[0].casefold()
            )
            enabled_family_gene_keys = render_matched_gene_table(
                matched_family_entries,
                len(family_keys),
                f"family_matched_gene_selection_{family_selection_id}",
                mean_expression_by_study(all_members_by_study, family_keys),
            )
        filtered_family_members = filter_matched_genes(
            all_members_by_study,
            enabled_family_gene_keys,
        )
        members_by_study = {
            key: filtered_family_members.get(key, members.iloc[0:0].copy())
            for key, members in all_members_by_study.items()
        }
        threshold = 1.0
        row_zscore = st.toggle(
            "Show relative pattern within each gene",
            value=True,
            help="Convert each gene's log₂(TPM + 1) values to z-scores across conditions. Positive means above that gene's average; negative means below it.",
        )

    if custom_family and family_queries:
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
    elif matched_family_entries and not enabled_family_gene_keys:
        st.info("Turn on at least one matched gene to show family results.")
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

            selected_names = ranking["gene"].tolist()
            selected_grouped = grouped_median(long, field)
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
                concise,
                hide_index=True,
                width="stretch",
                column_config={
                    "Mean TPM": st.column_config.NumberColumn(format="%.2f")
                },
            )
            with st.expander("Download family matrix"):
                download_tsv(
                    "Download complete family matrix",
                    matrix_for_genes(dataset, members),
                    f"{key}_{'custom' if custom_family else family_name.split('(')[-1].rstrip(')').casefold()}_family_tpm.tsv",
                    f"family_download_{key}",
                )

elif mode == "Differential expression":
    page_heading(
        "Differential expression",
        "Browse count-aware DESeq2 contrasts produced by nf-core/differentialabundance.",
    )
    st.caption(
        "The MA plot shows mean normalized count on the horizontal axis and target-versus-reference log₂ fold change on the vertical axis."
    )
    with st.container(key="differential_setup_panel"):
        comparison_key = st.selectbox(
            "Study",
            options=differential_study_keys,
            format_func=lambda key: (
                datasets[key].label
                if key in differential_contrasts
                else f"{datasets[key].label} — NOT AVAILABLE"
            ),
            key="differential_study",
        )
        comparison_dataset = datasets[comparison_key]
        available_contrasts = differential_contrasts.get(comparison_key, [])
        if available_contrasts:
            selected_contrast_id = st.selectbox(
                "Contrast · target vs reference",
                options=[contrast.contrast_id for contrast in available_contrasts],
                format_func=lambda contrast_id: contrast_label(
                    comparison_dataset,
                    next(
                        contrast
                        for contrast in available_contrasts
                        if contrast.contrast_id == contrast_id
                    ),
                ),
                key="differential_contrast",
            )
            selected_contrast = next(
                contrast
                for contrast in available_contrasts
                if contrast.contrast_id == selected_contrast_id
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

    if not available_contrasts:
        st.subheader("NOT AVAILABLE")
        st.info(
            "NOT AVAILABLE — this study has no bundled "
            "nf-core/differentialabundance results. The atlas does not compute "
            "differential-expression statistics from TPM values."
        )
    else:
        try:
            comparison_results = differential_results_cached(
                DATA_SCHEMA_VERSION,
                comparison_key,
                selected_contrast.contrast_id,
            )
        except ValueError as exc:
            st.warning(str(exc))
        else:
            target_label = condition_label(
                comparison_dataset,
                selected_contrast.variable,
                selected_contrast.target,
            )
            reference_label = condition_label(
                comparison_dataset,
                selected_contrast.variable,
                selected_contrast.reference,
            )
            target_samples, reference_samples = contrast_sample_counts(
                comparison_dataset, selected_contrast
            )
            plotted_results = comparison_results[comparison_results["ma_plot_eligible"]]
            significant_count = int((plotted_results["fdr"] < fdr_threshold).sum())
            omitted_count = len(comparison_results) - len(plotted_results)
            fold_change_range = differential_log2_fold_change_range(comparison_results)
            abundance_range = differential_abundance_range(comparison_results)
            off_scale_count = int(
                (
                    (plotted_results["log2_fold_change"] < fold_change_range[0])
                    | (plotted_results["log2_fold_change"] > fold_change_range[1])
                ).sum()
            )
            low_abundance_off_scale = int(
                (np.log10(plotted_results["base_mean"]) < abundance_range[0]).sum()
            )
            target_metric, reference_metric, significant_metric = st.columns(3)
            target_metric.metric("Target samples", target_samples, help=target_label)
            reference_metric.metric(
                "Reference samples", reference_samples, help=reference_label
            )
            significant_metric.metric(f"Colored genes · FDR < {fdr_threshold:g}", significant_count)
            st.plotly_chart(
                differential_ma_figure(comparison_results, fdr_threshold),
                width="stretch",
                key=f"differential_{comparison_key}_{selected_contrast.contrast_id}",
            )
            st.caption(
                f"{len(plotted_results):,} genes plotted. {omitted_count:,} genes with zero or unavailable DESeq2 base mean / fold change are omitted. The initial view excludes {low_abundance_off_scale:,} extreme low-abundance points and {off_scale_count:,} extreme fold changes; use Plotly zoom to inspect them."
            )
            st.caption(
                "Gray genes do not pass the selected FDR threshold. Significant genes are gold and drawn last, so gray points cannot cover them. All markers are fully opaque."
            )
            st.caption(
                f"Precomputed by nf-core/differentialabundance 2.0.0 with DESeq2 from Salmon gene-level counts. Positive log₂ fold change means higher in {target_label}; negative means higher in {reference_label}."
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
                    "base_mean": "DESeq2 base mean",
                    "log2_fold_change": "Log₂ fold change",
                    "fold_change": "Fold change (target / reference)",
                    "lfc_se": "Log₂ fold-change SE",
                    "p_value": "Raw p-value",
                    "fdr": "FDR",
                    "passes_fdr": threshold_column,
                }
            )[
                [
                    "Gene",
                    "Stable ID",
                    "DESeq2 base mean",
                    "Log₂ fold change",
                    "Fold change (target / reference)",
                    "Log₂ fold-change SE",
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
                    "DESeq2 base mean": st.column_config.NumberColumn(format="%.3f"),
                    "Log₂ fold change": st.column_config.NumberColumn(format="%.3f"),
                    "Fold change (target / reference)": st.column_config.NumberColumn(format="%.3f"),
                    "Log₂ fold-change SE": st.column_config.NumberColumn(format="%.3f"),
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
                    "base_mean": "DESeq2 base mean",
                    "log2_fold_change": "Log₂ fold change",
                    "fold_change": "Fold change (target / reference)",
                    "lfc_se": "Log₂ fold-change SE",
                    "p_value": "Raw p-value",
                    "fdr": "FDR",
                    "passes_fdr": threshold_column,
                }
            )[
                [
                    "Gene",
                    "Stable ID",
                    "DESeq2 base mean",
                    "Log₂ fold change",
                    "Fold change (target / reference)",
                    "Log₂ fold-change SE",
                    "Raw p-value",
                    "FDR",
                    threshold_column,
                ]
            ]
            download_table.insert(0, "Reference", reference_label)
            download_table.insert(0, "Target", target_label)
            download_table.insert(0, "Contrast", selected_contrast.contrast_id)
            download_table.insert(0, "Study", comparison_dataset.label)
            download_tsv(
                "Download all differential-expression results",
                download_table,
                f"{comparison_key}_{selected_contrast.contrast_id}_deseq2.tsv",
                f"differential_download_{comparison_key}_{selected_contrast.contrast_id}",
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
            available_cluster_genes = len(cluster_dataset.values)
            variable_genes = st.select_slider(
                "Genes used",
                options=cluster_gene_count_options(available_cluster_genes),
                value=available_cluster_genes,
                format_func=lambda value: (
                    f"All ({available_cluster_genes:,})"
                    if value == available_cluster_genes
                    else f"{value:,} most variable"
                ),
                help="All genes are used by default. Smaller choices keep the genes with the highest variance after log-transforming TPM.",
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
