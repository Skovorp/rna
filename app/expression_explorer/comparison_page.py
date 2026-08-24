"""Shared native Streamlit renderer for paper-vs-reprocessed TPM comparisons."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
import streamlit as st


# Matches app/.streamlit/config.toml.
BACKGROUND = "#0E1518"
TEXT = "#EDF5F2"
GRID = "#24323a"

# Streamlit rewrites Magma's exact near-black low stop in dark mode.
_NEAR_BLACK = {"#000004": "#010106", "#000005": "#010106", "#000003": "#010106"}


@dataclass(frozen=True)
class ComparisonPage:
    page_title: str
    heading: str
    introduction: str
    asset_dir: Path
    report_filename: str
    rebuild_command: str
    matched_genes_help: str
    data_caption: str
    pca_grouping: str


@st.cache_data(show_spinner=False)
def load_bundle(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _protect_low_stop(colorscale: object) -> object:
    if not isinstance(colorscale, (list, tuple)):
        return colorscale
    stops = []
    for stop in colorscale:
        if isinstance(stop, (list, tuple)) and len(stop) == 2:
            position, color = stop
            stops.append([position, _NEAR_BLACK.get(str(color).lower(), color)])
        else:
            return colorscale
    return stops


def figure(bundle: dict, name: str) -> go.Figure:
    """Rebuild one figure with the app's chrome and original data colors."""
    built = go.Figure(bundle["figures"][name])
    for trace in built.data:
        if trace.type == "heatmap":
            trace.colorscale = _protect_low_stop(trace.colorscale)
    built.update_layout(
        title_text="",
        paper_bgcolor=BACKGROUND,
        plot_bgcolor=BACKGROUND,
        font_color=TEXT,
        legend_font_color=TEXT,
        margin={"l": 70, "r": 40, "t": 50, "b": 60},
    )
    built.update_xaxes(gridcolor=GRID, zerolinecolor=GRID, linecolor=GRID)
    built.update_yaxes(gridcolor=GRID, zerolinecolor=GRID, linecolor=GRID)
    for annotation in built.layout.annotations or ():
        annotation.font.color = TEXT
    return built


def render_comparison_page(config: ComparisonPage) -> None:
    bundle_path = config.asset_dir / "figures.json"
    report_path = config.asset_dir / config.report_filename

    st.set_page_config(
        page_title=config.page_title,
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="collapsed",
    )
    st.markdown(
        """
        <style>
        .block-container p, .block-container li { font-size: 1.08rem; line-height: 1.65; }
        .block-container h2 { font-size: 1.7rem; margin-top: 2.2rem; }
        </style>
        """,
        unsafe_allow_html=True,
    )
    st.page_link("app.py", label="Back to expression explorer", icon="⬅️")
    st.title(config.heading)
    st.markdown(config.introduction)

    if not bundle_path.is_file():
        st.error(
            "The comparison bundle is missing. Regenerate it with "
            f"`{config.rebuild_command}`."
        )
        st.stop()

    bundle = load_bundle(str(bundle_path))
    summary = bundle["summary"]
    agreement = summary["agreement"]
    discordance = summary["discordance"]
    transitions = summary["zero_transitions"]
    identity = summary["sample_identity"]
    pca = summary["pca"]

    matched, correlation, error, within = st.columns(4)
    matched.metric(
        "Matched genes",
        f"{summary['matched_genes']:,}",
        help=config.matched_genes_help,
    )
    correlation.metric(
        "Pearson r",
        f"{agreement['pearson_log2_tpm_plus_1']:.3f}",
        help="On log₂(TPM + 1). Spearman is "
        f"{agreement['spearman_log2_tpm_plus_1']:.3f}.",
    )
    error.metric(
        "Median |log₂ error|",
        f"{agreement['median_absolute_log2_error']:.3f}",
        help="Median absolute difference in log₂(TPM + 1) across all matched "
        "gene-sample pairs.",
    )
    within.metric(
        "Within ±1 log₂",
        f"{agreement['fraction_abs_error_le_1']:.1%}",
        help="Share of gene-sample pairs agreeing to within one log₂ unit.",
    )
    st.caption(config.data_caption)

    analysis = st.segmented_control(
        "Analysis",
        (
            "TPM agreement",
            "Zero transitions",
            "Sample PCA",
            "Sample identity",
        ),
        default="TPM agreement",
        width="stretch",
        help="Only the selected analysis is loaded, so the page stays responsive.",
    )

    if analysis == "Zero transitions":
        st.markdown("## Exact zero and non-zero transitions")
        st.plotly_chart(figure(bundle, "zero_transition"), width="stretch", theme=None)
        st.markdown(
            f"**Green** marks published 0 becoming reprocessed non-zero; **red** marks "
            f"published non-zero becoming reprocessed 0. An exact zero can mean no "
            f"compatible fragments were assigned under that quantification model; "
            f"it is not a universal biological absence threshold. There are "
            f"{transitions['published_zero_to_reanalysis_nonzero_count']:,} pairs in "
            f"the first direction and "
            f"{transitions['published_nonzero_to_reanalysis_zero_count']:,} in the "
            f"second. At a non-zero-side threshold of 1 TPM these fall to "
            f"{transitions['published_zero_to_reanalysis_ge_1_count']:,} and "
            f"{transitions['published_ge_1_to_reanalysis_zero_count']:,}; at 10 TPM, "
            f"{transitions['published_zero_to_reanalysis_ge_10_count']:,} and "
            f"{transitions['published_ge_10_to_reanalysis_zero_count']:,}."
        )

        for table in bundle["zero_transition_tables"]:
            st.markdown(f"**{table['title']}**")
            st.dataframe(pd.DataFrame(table["rows"]), width="stretch", hide_index=True)
    elif analysis == "Sample PCA":
        st.markdown("## Sample PCA")
        st.plotly_chart(figure(bundle, "pca"), width="stretch", theme=None)
        st.markdown(
            f"PCA uses all {pca['genes_used']:,} one-to-one matched genes, with no "
            f"expression or variability cutoff. Expression is log-transformed and "
            f"each gene is standardized across samples. {config.pca_grouping} The "
            f"separate PCAs compare biological geometry; the joint PCA also exposes "
            f"method-specific shifts, where a circle is the published profile and a "
            f"cross is ours. Pairwise sample distances correlate at "
            f"{pca['pairwise_sample_distance_pearson']:.4f}."
        )
    elif analysis == "Sample identity":
        st.markdown("## Sample identity")
        st.plotly_chart(figure(bundle, "correlation"), width="stretch", theme=None)
        st.markdown(
            f"The diagonal compares the same biological sample across processing "
            f"methods. A diagonal maximum in each row argues against sample swaps, "
            f"and that holds for {identity['matching_sample_top1_count']} of "
            f"{summary['samples']} samples. Matching-sample correlations run from "
            f"{identity['min_matching_sample_correlation']:.4f} to "
            f"{identity['max_matching_sample_correlation']:.4f}."
        )
    else:
        st.markdown("## TPM agreement")
        st.plotly_chart(figure(bundle, "error"), width="stretch", theme=None)
        st.markdown(
            f"Errors are reprocessed minus published values. The diagonal bands are "
            f"forced by the coordinates: when published TPM is zero the error is "
            f"+2 × average log-expression, and when reprocessed TPM is zero it is "
            f"−2 × average log-expression. Of the "
            f"{discordance['abs_log2_error_gt_2_count']:,} pairs with absolute error "
            f"above 2, {discordance['severe_pairs_with_exact_zero_fraction']:.1%} "
            f"contain an exact zero. Density plots clip only the outer 0.1% for "
            f"readable axes; the summary metrics use the full distribution."
        )

    if report_path.is_file():
        st.divider()
        st.download_button(
            "Download the standalone report",
            data=report_path.read_text(encoding="utf-8"),
            file_name=config.report_filename,
            mime="text/html",
            help="A self-contained HTML copy of these figures for sharing. Its "
            "charts load plotly.js from a CDN, so it needs a network connection.",
        )
