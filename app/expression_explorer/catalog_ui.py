"""Responsive catalog rows with links to source reads and calculated TPM."""

from html import escape
from pathlib import Path

import streamlit as st
from streamlit import runtime

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
                if entry.pending_reprocessing:
                    explore.markdown("Not yet available")
                elif entry.comparison:
                    explore.markdown(f"[Published vs reprocessed]({entry.comparison})")
                elif entry.explorer:
                    # Same-tab links preserve a straightforward route into the app.
                    explore.markdown(
                        f'<a href="{entry.explorer}" target="_self">Open in Genes</a>',
                        unsafe_allow_html=True,
                    )
                else:
                    explore.markdown("No published counterpart")
                with downloads:
                    if entry.pending_reprocessing:
                        st.markdown("TPM processing pending")
                    elif entry.tpm_file:
                        path = tpm_download_path(entry.key, expression_dir, unlocked)
                        # Use the same session-managed file delivery as
                        # st.download_button, rendered as an ordinary link.
                        url = runtime.get_instance().media_file_mgr.add(
                            _file_bytes(str(path), path.stat().st_mtime_ns),
                            mimetype="application/gzip",
                            coordinates=f"catalog.tpm.{entry.key}",
                            file_name=path.name,
                            is_for_static_download=True,
                        )
                        st.markdown(
                            f'<a href="{escape(url, quote=True)}" '
                            f'download="{escape(path.name, quote=True)}" target="_self">'
                            'Download TPM tables</a>',
                            unsafe_allow_html=True,
                        )
                    else:
                        for source, url in entry.raw_sources:
                            st.markdown(f'[Download raw data]({url} "{source}")')
                            if len(entry.raw_sources) > 1:
                                st.caption(source)
        if not unlocked:
            with st.container(key="dataset_row_private_prompt"):
                st.button(
                    "Enter the password to unlock more private datasets",
                    on_click=on_unlock,
                    key="catalog_unlock_private",
                )
