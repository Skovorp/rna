"""Responsive catalog rows with native Streamlit download controls."""

from pathlib import Path

import streamlit as st

from .catalog import tpm_download_path, visible_catalog


@st.cache_resource(show_spinner=False, max_entries=8)
def _file_bytes(path: str, modified_ns: int) -> bytes:
    del modified_ns
    return Path(path).read_bytes()


def render_catalog(expression_dir: Path, unlocked: bool, on_unlock) -> None:
    st.markdown("## Datasets")
    st.markdown(
        """<style>
        .st-key-dataset_catalog [class*="st-key-dataset_row_"] {
            border-bottom: 1px solid rgba(148, 163, 184, .25);
            padding: .9rem 0;
        }
        .st-key-dataset_catalog [data-testid="stDownloadButton"] button {
            text-align: left;
        }
        @media (max-width: 640px) {
            .st-key-dataset_catalog_header { display: none; }
            .st-key-dataset_catalog [class*="st-key-dataset_row_"]
            [data-testid="stHorizontalBlock"] { gap: .35rem; }
        }
        </style>""",
        unsafe_allow_html=True,
    )
    widths = [2.4, 4.5, 2.1, 2.4]
    with st.container(key="dataset_catalog"):
        with st.container(key="dataset_catalog_header"):
            for column, label in zip(
                st.columns(widths),
                ("Dataset", "What is displayed", "Explore / compare", "Data downloads"),
            ):
                column.markdown(f"**{label}**")
        for entry in visible_catalog(unlocked):
            with st.container(key=f"dataset_row_{entry.key}"):
                name, description, explore, downloads = st.columns(widths)
                name.markdown(f"**{entry.label}**")
                description.markdown(entry.description)
                if entry.comparison:
                    explore.markdown(f"[Published vs reprocessed]({entry.comparison})")
                elif entry.explorer:
                    label = (
                        "Explore fruitless exons" if entry.key == "basrur"
                        else "Open UCSC atlas" if entry.key == "goldman"
                        else "Open in Genes"
                    )
                    # Same-tab links preserve a straightforward route into the app.
                    explore.markdown(
                        f'<a href="{entry.explorer}" target="_self">{label}</a>',
                        unsafe_allow_html=True,
                    )
                else:
                    explore.markdown("No published counterpart")
                with downloads:
                    for source, url in entry.raw_sources:
                        st.markdown(f'[Download raw data]({url} "{source}")')
                        if len(entry.raw_sources) > 1:
                            st.caption(source)
                    if not entry.raw_sources:
                        st.caption("Lab-provided raw reads")
                    if entry.tpm_file:
                        path = tpm_download_path(entry.key, expression_dir, unlocked)
                        st.download_button(
                            "Download TPM tables",
                            _file_bytes(str(path), path.stat().st_mtime_ns),
                            file_name=path.name,
                            mime="application/gzip",
                            key=f"download_tpm_{entry.key}",
                            help="Complete gene-by-sample TPM table (TSV, gzip compressed).",
                            on_click="ignore",
                        )
        if not unlocked:
            with st.container(key="dataset_row_private_prompt"):
                st.button(
                    "Enter the password to unlock more private datasets",
                    on_click=on_unlock,
                    key="catalog_unlock_private",
                )
