"""Load differential-expression results produced by nf-core pipelines."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from expression_explorer.data import ExpressionDataset


@dataclass(frozen=True)
class DifferentialContrast:
    contrast_id: str
    dataset_key: str
    variable: str
    reference: str
    target: str
    result_path: Path
    method: str


MANIFEST_COLUMNS = {
    "contrast_id",
    "dataset_key",
    "variable",
    "reference",
    "target",
    "result_file",
    "method",
}
RESULT_COLUMNS = (
    "gene_id",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "pvalue",
    "padj",
)


def load_differential_contrasts(
    expression_dir: Path | str,
) -> dict[str, list[DifferentialContrast]]:
    """Load the manifests for bundled nf-core differential results."""
    expression_dir = Path(expression_dir)
    manifests = sorted(expression_dir.glob("*_deseq2/contrasts.tsv"))
    contrasts_by_dataset: dict[str, list[DifferentialContrast]] = {}
    seen_ids: set[tuple[str, str]] = set()

    for manifest_path in manifests:
        manifest = pd.read_csv(manifest_path, sep="\t", dtype=str).fillna("")
        missing_columns = MANIFEST_COLUMNS - set(manifest.columns)
        if missing_columns:
            raise ValueError(
                f"Differential manifest is missing columns {sorted(missing_columns)}: "
                f"{manifest_path}"
            )
        for row in manifest.itertuples(index=False):
            identity = (row.dataset_key, row.contrast_id)
            if identity in seen_ids:
                raise ValueError(f"Duplicate differential contrast: {identity}")
            seen_ids.add(identity)
            result_path = manifest_path.parent / row.result_file
            if not result_path.is_file():
                raise FileNotFoundError(
                    f"Differential result is missing for {row.contrast_id}: {result_path}"
                )
            contrast = DifferentialContrast(
                contrast_id=row.contrast_id,
                dataset_key=row.dataset_key,
                variable=row.variable,
                reference=row.reference,
                target=row.target,
                result_path=result_path,
                method=row.method,
            )
            contrasts_by_dataset.setdefault(row.dataset_key, []).append(contrast)
    return contrasts_by_dataset


def condition_label(dataset: ExpressionDataset, variable: str, value: str) -> str:
    """Return the human-readable label associated with a contrast value."""
    samples = dataset.samples
    if variable not in samples.columns:
        raise ValueError(f"Unknown contrast variable: {variable}")
    matches = samples.loc[samples[variable].astype(str).eq(value)]
    if matches.empty:
        raise ValueError(f"Contrast value is absent from sample metadata: {value}")
    if "condition_label" in matches.columns:
        labels = [
            str(label).strip()
            for label in matches["condition_label"].drop_duplicates()
            if str(label).strip()
        ]
        if len(labels) == 1:
            return labels[0]
    return value


def contrast_label(dataset: ExpressionDataset, contrast: DifferentialContrast) -> str:
    """Format a target-versus-reference contrast for display."""
    target = condition_label(dataset, contrast.variable, contrast.target)
    reference = condition_label(dataset, contrast.variable, contrast.reference)
    return f"{target} vs {reference}"


def contrast_values(contrasts: list[DifferentialContrast]) -> list[str]:
    """Return every condition in stable manifest order."""
    return list(
        dict.fromkeys(
            value
            for contrast in contrasts
            for value in (contrast.reference, contrast.target)
        )
    )


def contrast_for_pair(
    contrasts: list[DifferentialContrast], target: str, reference: str
) -> DifferentialContrast:
    """Find the single precomputed result for an unordered condition pair."""
    if target == reference:
        raise ValueError("Target and reference must be different conditions")
    selected = [
        contrast
        for contrast in contrasts
        if {contrast.target, contrast.reference} == {target, reference}
    ]
    if len(selected) != 1:
        raise ValueError(
            f"Expected one differential contrast for {target} vs {reference}; "
            f"found {len(selected)}"
        )
    return selected[0]


def orient_differential_results(
    results: pd.DataFrame,
    contrast: DifferentialContrast,
    target: str,
    reference: str,
) -> pd.DataFrame:
    """Orient a stored pairwise result to the viewer's chosen direction.

    The bundle stores one DESeq2 result per unordered pair. When the selected
    target is the stored reference, reversing the sign of log2 fold change is
    exactly the requested target/reference orientation; inferential statistics,
    standard errors, and base means are unchanged.
    """
    if {target, reference} != {contrast.target, contrast.reference}:
        raise ValueError(
            f"{target} vs {reference} does not match {contrast.contrast_id}"
        )
    if target == contrast.target:
        return results
    oriented = results.copy()
    oriented["log2_fold_change"] = -oriented["log2_fold_change"]
    oriented["fold_change"] = np.exp2(oriented["log2_fold_change"])
    return oriented


def contrast_sample_counts(
    dataset: ExpressionDataset, contrast: DifferentialContrast
) -> tuple[int, int]:
    """Return target and reference sample counts from the dataset metadata."""
    if contrast.variable not in dataset.samples.columns:
        raise ValueError(f"Unknown contrast variable: {contrast.variable}")
    values = dataset.samples[contrast.variable].astype(str)
    return int(values.eq(contrast.target).sum()), int(values.eq(contrast.reference).sum())


def load_differential_results(
    dataset: ExpressionDataset,
    contrast: DifferentialContrast,
) -> pd.DataFrame:
    """Load one precomputed DESeq2 result and attach atlas gene names."""
    if contrast.dataset_key != dataset.key:
        raise ValueError(
            f"Contrast {contrast.contrast_id} belongs to {contrast.dataset_key}, "
            f"not {dataset.key}"
        )

    raw = pd.read_csv(contrast.result_path, sep="\t", compression="infer")
    missing_columns = set(RESULT_COLUMNS) - set(raw.columns)
    if missing_columns:
        raise ValueError(
            f"Differential result is missing columns {sorted(missing_columns)}: "
            f"{contrast.result_path}"
        )
    if raw["gene_id"].duplicated().any():
        raise ValueError(f"Differential result has duplicate gene IDs: {contrast.result_path}")

    annotations = (
        dataset.genes[["stable_id", "display_name"]]
        .drop_duplicates("stable_id")
        .set_index("stable_id")
    )
    results = raw[list(RESULT_COLUMNS)].copy()
    for column in ("baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"):
        results[column] = pd.to_numeric(results[column], errors="coerce")
    results["gene"] = results["gene_id"].map(annotations["display_name"])
    results["gene"] = results["gene"].fillna(results["gene_id"])
    results = results.rename(
        columns={
            "gene_id": "stable_id",
            "baseMean": "base_mean",
            "log2FoldChange": "log2_fold_change",
            "lfcSE": "lfc_se",
            "pvalue": "p_value",
            "padj": "fdr",
        }
    )
    results["fold_change"] = np.exp2(results["log2_fold_change"])
    results["ma_plot_eligible"] = (
        results["base_mean"].gt(0)
        & np.isfinite(results["base_mean"])
        & np.isfinite(results["log2_fold_change"])
    )
    results = results[
        [
            "gene",
            "stable_id",
            "base_mean",
            "log2_fold_change",
            "fold_change",
            "lfc_se",
            "p_value",
            "fdr",
            "ma_plot_eligible",
        ]
    ]
    return results.sort_values(
        ["fdr", "p_value", "log2_fold_change"],
        ascending=[True, True, False],
        na_position="last",
        kind="stable",
    ).reset_index(drop=True)
