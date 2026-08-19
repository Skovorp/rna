from __future__ import annotations

from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components


# One report, not two: the full report already embeds the same zero-transition
# figure as the old standalone page, plus the per-gene transition tables that
# page lacked. The generator still writes the standalone file; the app ignores
# it.
#
# Embedded with components.html rather than served from app/static/, because
# Streamlit's static server returns .html as text/plain with nosniff and an
# iframe pointed at it would show raw source. Embedded size therefore matters:
# scripts/theme_comparison_reports.py strips the inlined ~4.9 MB plotly.js to a
# CDN reference, taking this file well under 0.5 MB.
ASSET_DIR = Path(__file__).resolve().parents[1] / "assets" / "ovary_comparison"
REPORT = "elife_ovary_tpm_full_report.html"

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

st.markdown(
    "The report below covers, in order: **overall TPM agreement** and the "
    "error distribution, **exact zero ↔ non-zero transitions** with the "
    "genes driving them, **PCA** of the two matrices, and a **sample-identity "
    "correlation** check against sample swaps."
)


@st.cache_data(show_spinner=False)
def report_html(file_name: str) -> str:
    return (ASSET_DIR / file_name).read_text(encoding="utf-8")


if not (ASSET_DIR / REPORT).is_file():
    st.error(f"Bundled report is missing: {REPORT}")
else:
    html = report_html(REPORT)
    components.html(html, height=3200, scrolling=True)
    st.download_button(
        "Download this report",
        data=html,
        file_name=REPORT,
        mime="text/html",
        help="Charts in the downloaded file load plotly.js from a CDN, so it "
        "needs a network connection to render.",
    )
