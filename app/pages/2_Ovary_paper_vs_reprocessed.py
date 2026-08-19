from __future__ import annotations

from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components


# Reports are embedded with components.html rather than served from app/static/:
# Streamlit's static server returns .html as text/plain with nosniff, so an
# iframe pointed at it shows raw source. That makes the embedded size the thing
# that matters — scripts/theme_comparison_reports.py strips the inlined ~4.9 MB
# plotly.js down to a CDN reference, taking these files to well under 0.5 MB.
ASSET_DIR = Path(__file__).resolve().parents[1] / "assets" / "ovary_comparison"

st.set_page_config(
    page_title="Ovary paper vs reprocessed · Aedes RNA Atlas",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed",
)

st.page_link("app.py", label="Back to expression explorer", icon="⬅️")
st.title("Ovary · paper vs reprocessed")
st.markdown(
    "Comparison of the **published TPM supplement** from Venkataraman et al., "
    "eLife 2023 against **our reprocessing** of the same raw reads "
    "(`PRJNA796320`). Both matrices cover the same 33 ovary samples; genes are "
    "matched through an identifier crosswalk. See the **Methods** page for the "
    "exact pipeline and parameters used for every reprocessed dataset."
)

agreement, discordance, transitions = st.columns(3)
agreement.metric(
    "Correlation",
    "0.972",
    help="Pearson r on log₂(TPM + 1) across all 604,032 matched gene-sample "
    "pairs. Spearman is 0.936.",
)
discordance.metric(
    "Pairs disagreeing > 2 log₂",
    "2.8%",
    help="17,007 of 604,032 gene-sample pairs. 1,543 genes disagree that much "
    "in at least one sample; 112 do so in every sample.",
)
transitions.metric(
    "Published ≥ 10 TPM, ours 0",
    "4,073",
    help="The asymmetry runs one way: only 163 pairs are ≥ 10 TPM in our "
    "matrix while exactly zero in the published one.",
)
st.caption(
    "Regenerated from the current reprocessed matrix, so these figures track "
    "what the atlas actually displays."
)

REPORTS = {
    "TPM agreement report": (
        "elife_ovary_tpm_full_report.html",
        2400,
        "Per-sample and per-gene agreement between published and reprocessed "
        "TPM values: correlations, log2 error distributions, PCA, and the "
        "most discordant genes.",
    ),
    "Zero ↔ non-zero transitions": (
        "elife_ovary_zero_nonzero_transitions.html",
        1400,
        "Genes that are exactly zero in one matrix but expressed in the "
        "other — the dominant driver of severe disagreements.",
    ),
}


@st.cache_data(show_spinner=False)
def report_html(file_name: str) -> str:
    return (ASSET_DIR / file_name).read_text(encoding="utf-8")


choice = st.radio("Report", list(REPORTS), horizontal=True)
file_name, height, description = REPORTS[choice]
st.caption(description)

if not (ASSET_DIR / file_name).is_file():
    st.error(f"Bundled report is missing: {file_name}")
else:
    html = report_html(file_name)
    components.html(html, height=height, scrolling=True)
    st.download_button(
        "Download this report",
        data=html,
        file_name=file_name,
        mime="text/html",
        help="Charts in the downloaded file load plotly.js from a CDN, so it "
        "needs a network connection to render.",
    )
