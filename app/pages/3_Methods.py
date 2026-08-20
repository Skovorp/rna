from __future__ import annotations

from pathlib import Path

import streamlit as st


METHODS_PATH = (
    Path(__file__).resolve().parents[2] / "expression" / "METHODS.md"
)

st.set_page_config(
    page_title="Methods - Aedes RNA Atlas",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed",
)

st.markdown(
    """
    <style>
    .block-container { max-width: 1100px; }
    /* Long-form reference page: read at a comfortable size, not chart-label size. */
    .block-container p,
    .block-container li,
    .block-container td,
    .block-container th { font-size: 1.15rem; line-height: 1.7; }
    .block-container h2 { font-size: 1.9rem; margin-top: 2rem; }
    .block-container h3 { font-size: 1.45rem; }
    .block-container code { font-size: 1.02rem; }
    .block-container pre code { font-size: 0.98rem; line-height: 1.55; }
    </style>
    """,
    unsafe_allow_html=True,
)

st.page_link("app.py", label="Back to expression explorer", icon="⬅️")
st.title("Reprocessing methods")
st.markdown(
    "Pipeline, reference, and parameters used for every dataset we "
    "reprocessed from raw reads (ovary, midgut, crop). Paper datasets shown "
    "in the atlas keep their published values and are not covered by these "
    "methods."
)

if not METHODS_PATH.is_file():
    st.error("METHODS.md is missing from the expression bundle.")
else:
    st.markdown(METHODS_PATH.read_text(encoding="utf-8"))
