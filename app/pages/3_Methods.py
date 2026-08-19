from __future__ import annotations

from pathlib import Path

import streamlit as st


METHODS_PATH = (
    Path(__file__).resolve().parents[2] / "expression" / "METHODS.md"
)

st.set_page_config(
    page_title="Methods · Aedes RNA Atlas",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed",
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
