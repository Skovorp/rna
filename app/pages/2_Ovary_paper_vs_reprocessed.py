from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
import streamlit as st


# Rendered as ordinary page content from a plotly JSON bundle, not as an
# embedded HTML report. Embedding meant an iframe with its own inner scrollbar,
# its own copy of plotly, and the whole file crossing the websocket on every
# rerun. Regenerate the bundle with scripts/rebuild_ovary_comparison.py.
ASSET_DIR = Path(__file__).resolve().parents[1] / "assets" / "ovary_comparison"
BUNDLE = ASSET_DIR / "figures.json"
REPORT = "elife_ovary_tpm_full_report.html"

st.set_page_config(
    page_title="Ovary paper vs reprocessed - Aedes RNA Atlas",
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


@st.cache_data(show_spinner=False)
def load_bundle() -> dict:
    return json.loads(BUNDLE.read_text(encoding="utf-8"))


# Matches app/.streamlit/config.toml.
BACKGROUND = "#0E1518"
TEXT = "#EDF5F2"
GRID = "#24323a"

# Substituted for Magma's near-black low stop. Streamlit rewrites chart colours
# for contrast in dark mode by exact value, and it maps #000004 to a bright
# purple, which washed the density panel out. This is visually identical black
# but not the value it looks for.
_NEAR_BLACK = {"#000004": "#010106", "#000005": "#010106", "#000003": "#010106"}


def _protect_low_stop(colorscale: object) -> object:
    """Nudge near-black colorscale stops off the value Streamlit rewrites."""
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
    """Rebuild one figure, themed to the app and stripped of its report title.

    title_text="" rather than title=None: the latter renders the string
    "undefined" above the chart. Only the chrome is recolored, never the
    template: switching to plotly_dark replaces the default sequential
    colorscale and turns the density heatmaps bright magenta.

    Heatmap colorscales are restated by name because an expanded stop list gets
    partially rewritten before it reaches the browser: Magma's black low end
    arrived as #ab63fa, washing the density panel out to bright purple.
    """
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


st.page_link("app.py", label="Back to expression explorer", icon="⬅️")
st.title("Ovary: paper vs reprocessed")
st.markdown(
    "Comparison of the **published TPM supplement** from Venkataraman et al., "
    "eLife 2023 against **our reprocessing** of the same raw reads "
    "(`PRJNA796320`). Both matrices cover the same 33 ovary samples; genes are "
    "matched through an identifier crosswalk. See the [Methods page](/Methods) "
    "for the exact pipeline and parameters used for every reprocessed dataset."
)

if not BUNDLE.is_file():
    st.error(
        "The comparison bundle is missing. Regenerate it with "
        "scripts/rebuild_ovary_comparison.py."
    )
    st.stop()

bundle = load_bundle()
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
    help=f"Out of {summary['published_genes']:,} genes in the published matrix.",
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

st.caption(
    "Regenerated from the current reprocessed matrix, so these figures track "
    "what the atlas actually displays."
)

st.markdown("## TPM agreement")
st.plotly_chart(figure(bundle, "error"), use_container_width=True, theme=None)
st.markdown(
    f"Errors are reprocessed minus published values. The diagonal bands are "
    f"forced by the coordinates: when published TPM is zero the error is "
    f"+2 x average log-expression, and when reprocessed TPM is zero it is "
    f"-2 x average log-expression. Of the "
    f"{discordance['abs_log2_error_gt_2_count']:,} pairs with absolute error "
    f"above 2, {discordance['severe_pairs_with_exact_zero_fraction']:.1%} "
    f"contain an exact zero. Density plots clip only the outer 0.1% for "
    f"readable axes; the summary metrics use the full distribution."
)

st.markdown("## Exact zero and non-zero transitions")
st.plotly_chart(
    figure(bundle, "zero_transition"), use_container_width=True, theme=None
)
st.markdown(
    f"**Green** marks published 0 becoming reprocessed non-zero; **red** marks "
    f"published non-zero becoming reprocessed 0. An exact zero can mean no "
    f"compatible fragments were assigned under that quantification model; it is "
    f"not a universal biological absence threshold. There are "
    f"{transitions['published_zero_to_reanalysis_nonzero_count']:,} pairs in the "
    f"first direction and "
    f"{transitions['published_nonzero_to_reanalysis_zero_count']:,} in the "
    f"second. At a non-zero-side threshold of 1 TPM these fall to "
    f"{transitions['published_zero_to_reanalysis_ge_1_count']:,} and "
    f"{transitions['published_ge_1_to_reanalysis_zero_count']:,}; at 10 TPM, "
    f"{transitions['published_zero_to_reanalysis_ge_10_count']:,} and "
    f"{transitions['published_ge_10_to_reanalysis_zero_count']:,}."
)

for table in bundle["zero_transition_tables"]:
    st.markdown(f"**{table['title']}**")
    st.dataframe(
        pd.DataFrame(table["rows"]), use_container_width=True, hide_index=True
    )

st.markdown("## Sample PCA")
st.plotly_chart(figure(bundle, "pca"), use_container_width=True, theme=None)
st.markdown(
    f"PCA uses all {pca['genes_used']:,} one-to-one matched genes, with no "
    f"expression or variability cutoff. Expression is log-transformed and each "
    f"gene is standardized across samples, matching the atlas PCA convention. "
    f"The separate PCAs compare biological geometry; the joint PCA also exposes "
    f"method-specific shifts, where a circle is the published profile and a "
    f"cross is ours. Pairwise sample distances correlate at "
    f"{pca['pairwise_sample_distance_pearson']:.4f}."
)

st.markdown("## Sample identity")
st.plotly_chart(
    figure(bundle, "correlation"), use_container_width=True, theme=None
)
st.markdown(
    f"The diagonal compares the same biological sample across processing "
    f"methods. A diagonal maximum in each row argues against sample swaps, and "
    f"that holds for "
    f"{identity['matching_sample_top1_count']} of {summary['samples']} samples. "
    f"Matching-sample correlations run from "
    f"{identity['min_matching_sample_correlation']:.4f} to "
    f"{identity['max_matching_sample_correlation']:.4f}."
)

if (ASSET_DIR / REPORT).is_file():
    st.divider()
    st.download_button(
        "Download the standalone report",
        data=(ASSET_DIR / REPORT).read_text(encoding="utf-8"),
        file_name=REPORT,
        mime="text/html",
        help="A self-contained HTML copy of these figures for sharing. Its "
        "charts load plotly.js from a CDN, so it needs a network connection.",
    )
