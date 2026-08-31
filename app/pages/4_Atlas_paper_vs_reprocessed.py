from pathlib import Path

from expression_explorer.comparison_page import ComparisonPage, render_comparison_page


render_comparison_page(
    ComparisonPage(
        page_title="Neurotranscriptome published vs reprocessed - Aedes RNA Atlas",
        heading="Neurotranscriptome: published vs reprocessed",
        introduction=(
            "Comparison of the **published AaegL.RU TPM matrix** from Matthews et al., "
            "BMC Genomics 2016 against **our reprocessing** of the same raw reads "
            "(`PRJNA236239`). The comparison uses all 122 samples present in the "
            "published matrix. Three additional recovered libraries shown in our "
            "reprocessed neurotranscriptome matrix are excluded here because they "
            "have no published TPM "
            "counterpart. Repeated historical gene identifiers in the paper matrix "
            "are collapsed by summing TPM before one-to-one identifier matching. See "
            "the [Methods page](/Methods) for the shared reprocessing pipeline."
        ),
        asset_dir=Path(__file__).resolve().parents[1] / "assets" / "atlas_comparison",
        report_filename="matthews_2016_atlas_tpm_full_report.html",
        rebuild_command="scripts/rebuild_atlas_comparison.py",
        matched_genes_help=(
            "Direct one-to-one AaegL.RU identifiers shared by the published and "
            "reprocessed gene-level matrices."
        ),
        data_caption=(
            "The published matrix has 16,154 unique historical gene identifiers; "
            "the reprocessed matrix uses the current 19,920-gene reference."
        ),
        pca_grouping=(
            "Points are colored by tissue; sample names retain feeding and "
            "reproductive-state detail."
        ),
    )
)
