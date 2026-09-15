"""Basrur's published exon measurements, kept in normalized-count units."""

from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st


@st.cache_data(show_spinner=False)
def load_fruitless(path: str, modified_ns: int) -> pd.DataFrame:
    del modified_ns
    return pd.read_csv(path, sep="\t")


def render_fruitless(expression_dir: Path) -> None:
    st.markdown("# Fruitless exon explorer")
    st.markdown(
        "Published *Aedes aegypti* **fruitless (fru) exon counts** from "
        "[Basrur et al. (2020)](https://elifesciences.org/articles/63982), "
        "Figure 1G–H, using the Matthews et al. RNA-seq samples. "
        "Values are normalized exon counts, not gene TPM."
    )
    path = expression_dir / "basrur_2020_fruitless_exon_counts.tsv.gz"
    data = load_fruitless(str(path), path.stat().st_mtime_ns)
    panel_label = st.radio(
        "Published expression table",
        ("Across tissues (m, f, c1 exons)", "Brain (all exons)"),
        horizontal=True,
    )
    panel = "Figure 1H" if panel_label.startswith("Across") else "Figure 1G"
    selected = data[data.panel.eq(panel)].copy()
    sex = st.multiselect("Sex", ["female", "male"], default=["female", "male"], key="fruitless_sex")
    selected = selected[selected.sex.isin(sex)]
    if panel == "Figure 1H":
        exons = st.multiselect("Exons", ["m", "f", "c1"], default=["m", "f", "c1"], key="fruitless_exons")
        selected = selected[selected.exon.isin(exons)]
        x, facet = "tissue", "exon"
    else:
        selected = selected.sort_values(["exon_number", "sex", "replicate"])
        x, facet = "exon", None
    if selected.empty:
        st.info("Choose at least one sex and exon to display the published measurements.")
    else:
        chart = px.strip(
            selected, x=x, y="normalized_count", color="sex", facet_row=facet,
            hover_data=["tissue", "exon", "replicate"], stripmode="group",
            category_orders={"sex": ["female", "male"], "exon": selected.exon.drop_duplicates().tolist()},
            labels={"normalized_count": "Normalized exon count", "tissue": "Tissue", "exon": "Exon", "sex": "Sex"},
            color_discrete_map={"female": "#ec4899", "male": "#53d6a5"},
        )
        chart.update_layout(height=300 * selected.exon.nunique() if facet else 480)
        chart.update_xaxes(tickangle=-45)
        st.plotly_chart(chart, use_container_width=True)
        st.caption("Each point is one published replicate. Missing replicates remain absent; they are not zeros.")
        st.dataframe(selected.drop(columns="exon_number"), hide_index=True, use_container_width=True)
    st.download_button(
        "Download published exon tables", data.to_csv(sep="\t", index=False).encode(),
        file_name="basrur_2020_fruitless_exon_counts.tsv", mime="text/tab-separated-values",
        on_click="ignore",
    )
    st.markdown("[Download the original Figure 1 workbook](https://cdn.elifesciences.org/articles/63982/elife-63982-fig1-data1-v3.xlsx)")
