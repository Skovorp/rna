from pathlib import Path

from expression_explorer.comparison_page import ComparisonPage, render_comparison_page


render_comparison_page(
    ComparisonPage(
        page_title="Ovary published vs reprocessed - Aedes RNA Atlas",
        heading="Ovary: published vs reprocessed",
        introduction=(
            "Comparison of the **published TPM supplement** from Venkataraman et al., "
            "eLife 2023 against **our reprocessing** of the same raw reads "
            "(`PRJNA796320`). Both matrices cover the same 33 ovary samples; genes "
            "are matched through an identifier crosswalk. See the [Methods page](/Methods) "
            "for the exact pipeline and parameters used for every reprocessed dataset."
        ),
        asset_dir=Path(__file__).resolve().parents[1] / "assets" / "ovary_comparison",
        report_filename="elife_ovary_tpm_full_report.html",
        rebuild_command="scripts/rebuild_ovary_comparison.py",
        matched_genes_help="Out of 18,304 genes in the published matrix.",
        data_caption=(
            "Regenerated from the current reprocessed matrix, so these figures track "
            "what the atlas actually displays."
        ),
        pca_grouping="Points are colored by reproductive state.",
    )
)
