from __future__ import annotations

from pathlib import Path

import streamlit as st
from expression_explorer.data import _cache_signature
from expression_explorer.gene_aliases import PublishedAliases
from expression_explorer.mapping_audit import mapping_archive, mapping_source_summary


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
    /* Keep the reference readable, with a compact source breakdown. */
    .block-container p,
    .block-container li { font-size: 1.15rem; line-height: 1.7; }
    .block-container td,
    .block-container th { font-size: 1.02rem; line-height: 1.5; }
    .block-container h2 { font-size: 1.9rem; margin-top: 2rem; }
    .block-container h3 { font-size: 1.45rem; }
    .block-container code { font-size: 1.02rem; }
    .block-container pre code { font-size: 0.98rem; line-height: 1.55; }
    </style>
    """,
    unsafe_allow_html=True,
)

st.page_link("app.py", label="Back to expression explorer", icon="⬅️")
st.title("Methods and gene mappings")

if not METHODS_PATH.is_file():
    st.error("METHODS.md is missing from the expression bundle.")
else:
    st.markdown(METHODS_PATH.read_text(encoding="utf-8"))

    registry = PublishedAliases.load(METHODS_PATH.parent)
    counts = mapping_source_summary(registry, METHODS_PATH.parent).set_index("source")
    goldman = counts.loc["goldman_2025_s1.4"]
    matthews = counts.loc["matthews_2016_s5"]
    morita = counts.loc["morita_2025"]
    jove = counts.loc["jove_2020"]
    venkataraman = counts.loc["venkataraman_2023"]

    def citation(key: str, label: str) -> str:
        return f"[{label}]({registry.sources[key]['citation']})"

    st.markdown("### Where the mappings come from")
    st.markdown(
        "| Source | Links supplied |\n| --- | --- |\n"
        f"| {citation('goldman_2025_s1.4', 'Goldman S1.4')} | "
        f"**{goldman.loc_aael_pairs:,} LOC ↔ AAEL pairs**; "
        f"{goldman.other_identifier_pairs:,} other ID pairs (e.g. Orco ↔ AAEL). "
        f"{goldman.distinct_names:,} distinct names across {goldman.source_rows:,} atlas rows. |\n"
        f"| {citation('matthews_2016_s5', 'Matthews 2016, S5')} | "
        f"**{matthews.id_name_pairs:,} ID ↔ name pairs** across {matthews.primary_identifiers:,} IDs "
        "(AAEL, RU and two combined AAEL labels). |\n"
        f"| {citation('morita_2025', 'Morita 2025')} | "
        f"**{morita.id_name_pairs:,} AAEL ↔ name pairs**; "
        f"{morita.id_name_pairs_also_in_goldman:,} also in Goldman, "
        f"{morita.id_name_pairs - morita.id_name_pairs_only_in_this_source - morita.id_name_pairs_also_in_goldman:,} "
        f"more in Matthews, {morita.id_name_pairs_only_in_this_source:,} only in Morita. |\n"
        f"| {citation('jove_2020', 'Jové 2020')} | "
        f"**{jove.dataset_local_gene_loc_pairs:,} local gene… ↔ LOC pairs** "
        f"(this dataset only); {jove.id_name_pairs:,} receptor IDs ↔ full names. |\n"
        f"| {citation('venkataraman_2023', 'Venkataraman 2023')} | "
        f"**{venkataraman.id_name_pairs:,} additional name pairs**: IDs and Symbols are identical; "
        "original IDs remain searchable. |"
    )
    st.caption("Pairs are deduplicated within each source. Sources overlap; these counts are not unique genes or dataset coverage.")
    st.markdown(
        f"**PPK names:** Goldman supplies {goldman.ppk_names:,} modern aliases (e.g. `ppk317`); "
        f"Matthews supplies {matthews.ppk_names:,} historical aliases (e.g. `ppk00873`); "
        f"Morita supplies {morita.ppk_names:,}, largely overlapping Goldman. "
        "For this gene, the published links give **`ppk00873 ↔ AAEL000873 ↔ LOC5567199 ↔ ppk317`**. "
        "Search works from either end; ambiguous names return separate genes."
    )

    @st.cache_data(show_spinner="Preparing the full gene mapping…", max_entries=4)
    def mapping_download(signature: str, unlocked: bool) -> bytes:
        return mapping_archive(METHODS_PATH.parent, unlocked)

    st.markdown("## Download the full gene mapping")
    st.write(
        "Every gene, including unmatched rows: IDs, aliases, both descriptions and source evidence. "
        "Also includes source counts, comparison pairs and full processing methods."
    )
    unlocked = bool(st.session_state.get("private_datasets_unlocked", False))
    signature = _cache_signature(METHODS_PATH.parent)
    for path in (METHODS_PATH.parent / "PROCESSING_DETAILS.md",
                 Path(__file__).resolve().parents[1] / "expression_explorer" / "mapping_audit.py"):
        signature += f":{path.stat().st_mtime_ns}"
    # Comparison updates invalidate the download independently of matrix inputs.
    for path in sorted((METHODS_PATH.parents[1] / "analysis" / "results").glob("*/matched_genes.tsv.gz")):
        signature += f":{path.stat().st_mtime_ns}"
    st.download_button(
        "Download full gene mapping (ZIP)",
        mapping_download(signature, unlocked),
        file_name="aedes_gene_mapping.zip", mime="application/zip",
    )
    st.caption("Private study rows are included when private datasets are unlocked." if unlocked
               else "Includes all public datasets. Unlock private datasets to include their rows.")
