from __future__ import annotations

from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components


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

REPORTS = {
    "TPM agreement report": (
        "elife_ovary_tpm_full_report.html",
        "Per-sample and per-gene agreement between published and reprocessed "
        "TPM values: correlations, log2 error distributions, PCA, and the "
        "most discordant genes.",
    ),
    "Zero ↔ non-zero transitions": (
        "elife_ovary_zero_nonzero_transitions.html",
        "Genes that are exactly zero in one matrix but expressed in the "
        "other — the dominant driver of severe disagreements.",
    ),
}

choice = st.radio("Report", list(REPORTS), horizontal=True)
file_name, description = REPORTS[choice]
st.caption(description)

report_path = ASSET_DIR / file_name
if not report_path.is_file():
    st.error(f"Bundled report is missing: {report_path.name}")
else:
    html = report_path.read_text(encoding="utf-8")
    components.html(html, height=2400, scrolling=True)
    st.download_button(
        "Download this report (self-contained HTML)",
        data=html,
        file_name=file_name,
        mime="text/html",
    )
