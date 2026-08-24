#!/usr/bin/env python3
"""Compare published and independently reanalysed eLife ovary TPM matrices."""

from __future__ import annotations

import argparse
import base64
from collections import Counter, defaultdict
from html import escape
import json
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
from scipy.spatial import procrustes
from scipy.spatial.distance import pdist
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import PCA


SOURCE_COLORS = {"Published": "#2563eb", "Reanalysis": "#dc2626"}
CONDITION_COLORS = [
    "#2563eb",
    "#f97316",
    "#16a34a",
    "#dc2626",
    "#7c3aed",
    "#0891b2",
    "#ca8a04",
    "#db2777",
    "#4f46e5",
    "#059669",
    "#9333ea",
]


def load_matrix(path: Path) -> tuple[pd.DataFrame, list[str]]:
    frame = pd.read_csv(path, sep="\t", dtype={"IDs": str, "Symbols": str})
    if list(frame.columns[:2]) != ["IDs", "Symbols"]:
        raise ValueError(f"{path} must begin with IDs and Symbols columns")
    if frame["IDs"].duplicated().any():
        raise ValueError(f"{path} contains duplicate IDs")
    samples = list(frame.columns[2:])
    numeric = frame[samples].apply(pd.to_numeric, errors="raise")
    if not np.isfinite(numeric.to_numpy()).all() or (numeric.to_numpy() < 0).any():
        raise ValueError(f"{path} contains invalid TPM values")
    frame[samples] = numeric
    return frame, samples


def build_gene_map(
    published: pd.DataFrame,
    reanalysis: pd.DataFrame,
    crosswalk: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, int]]:
    target_ids = set(reanalysis["IDs"])
    aliases: dict[str, set[str]] = defaultdict(set)
    for row in crosswalk.itertuples(index=False):
        if row.target_gene_id not in target_ids:
            continue
        for alias in (row.source_gene_id, row.source_gene_name):
            alias = str(alias).strip()
            if alias:
                aliases[alias].add(row.target_gene_id)

    candidates_by_source: dict[str, set[str]] = {}
    method_by_source: dict[str, str] = {}
    ambiguous_aliases = 0
    for source_id in published["IDs"]:
        candidates = set(aliases.get(source_id, set()))
        if source_id in target_ids:
            candidates.add(source_id)
        if len(candidates) == 1:
            candidates_by_source[source_id] = candidates
            method_by_source[source_id] = (
                "direct_identifier"
                if candidates == {source_id}
                else "exact_coordinate_crosswalk"
            )
        elif len(candidates) > 1:
            ambiguous_aliases += 1

    target_use = Counter(next(iter(value)) for value in candidates_by_source.values())
    rows = []
    duplicate_targets = 0
    published_lookup = published.set_index("IDs")
    reanalysis_lookup = reanalysis.set_index("IDs")
    for source_id in published["IDs"]:
        candidates = candidates_by_source.get(source_id)
        if not candidates:
            continue
        target_id = next(iter(candidates))
        if target_use[target_id] != 1:
            duplicate_targets += 1
            continue
        rows.append(
            {
                "published_id": source_id,
                "published_symbol": published_lookup.at[source_id, "Symbols"],
                "reanalysis_id": target_id,
                "reanalysis_symbol": reanalysis_lookup.at[target_id, "Symbols"],
                "mapping_method": method_by_source[source_id],
            }
        )
    gene_map = pd.DataFrame(rows)
    if gene_map.empty:
        raise ValueError("The annotation crosswalk produced no comparable genes")
    diagnostics = {
        "ambiguous_source_aliases": ambiguous_aliases,
        "sources_removed_for_duplicate_target": duplicate_targets,
    }
    return gene_map, diagnostics


def safe_pearson(x: np.ndarray, y: np.ndarray) -> float:
    if np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(pearsonr(x, y).statistic)


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    if np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(spearmanr(x, y).statistic)


def error_metrics(
    published_log: np.ndarray,
    reanalysis_log: np.ndarray,
    samples: list[str],
) -> tuple[pd.DataFrame, dict[str, float]]:
    residual = reanalysis_log - published_log
    rows = []
    for index, sample in enumerate(samples):
        sample_residual = residual[:, index]
        rows.append(
            {
                "sample": sample,
                "pearson_log2_tpm_plus_1": safe_pearson(
                    published_log[:, index], reanalysis_log[:, index]
                ),
                "spearman_log2_tpm_plus_1": safe_spearman(
                    published_log[:, index], reanalysis_log[:, index]
                ),
                "median_signed_log2_error": np.median(sample_residual),
                "median_absolute_log2_error": np.median(np.abs(sample_residual)),
                "rmse_log2_tpm_plus_1": np.sqrt(np.mean(sample_residual**2)),
                "error_q10": np.quantile(sample_residual, 0.10),
                "error_q25": np.quantile(sample_residual, 0.25),
                "error_q75": np.quantile(sample_residual, 0.75),
                "error_q90": np.quantile(sample_residual, 0.90),
                "fraction_abs_error_le_0_5": np.mean(np.abs(sample_residual) <= 0.5),
                "fraction_abs_error_le_1": np.mean(np.abs(sample_residual) <= 1.0),
            }
        )
    per_sample = pd.DataFrame(rows)
    flat_residual = residual.ravel()
    overall = {
        "pearson_log2_tpm_plus_1": safe_pearson(
            published_log.ravel(), reanalysis_log.ravel()
        ),
        "spearman_log2_tpm_plus_1": safe_spearman(
            published_log.ravel(), reanalysis_log.ravel()
        ),
        "median_signed_log2_error": float(np.median(flat_residual)),
        "median_absolute_log2_error": float(np.median(np.abs(flat_residual))),
        "rmse_log2_tpm_plus_1": float(np.sqrt(np.mean(flat_residual**2))),
        "error_q01": float(np.quantile(flat_residual, 0.01)),
        "error_q10": float(np.quantile(flat_residual, 0.10)),
        "error_q25": float(np.quantile(flat_residual, 0.25)),
        "error_q75": float(np.quantile(flat_residual, 0.75)),
        "error_q90": float(np.quantile(flat_residual, 0.90)),
        "error_q99": float(np.quantile(flat_residual, 0.99)),
        "fraction_abs_error_le_0_5": float(np.mean(np.abs(flat_residual) <= 0.5)),
        "fraction_abs_error_le_1": float(np.mean(np.abs(flat_residual) <= 1.0)),
    }
    return per_sample, overall


def discordance_metrics(
    published_tpm: np.ndarray,
    reanalysis_tpm: np.ndarray,
    published_log: np.ndarray,
    reanalysis_log: np.ndarray,
) -> dict[str, float | int]:
    absolute_error = np.abs(reanalysis_log - published_log)
    severe = absolute_error > 2.0
    exact_zero = (published_tpm == 0) | (reanalysis_tpm == 0)
    published_high_reanalysis_zero = (published_tpm >= 10) & (reanalysis_tpm == 0)
    reanalysis_high_published_zero = (reanalysis_tpm >= 10) & (published_tpm == 0)
    severe_samples_per_gene = severe.sum(axis=1)
    return {
        "gene_sample_pairs": int(absolute_error.size),
        "abs_log2_error_gt_2_count": int(severe.sum()),
        "abs_log2_error_gt_2_fraction": float(severe.mean()),
        "severe_pairs_with_exact_zero_count": int((severe & exact_zero).sum()),
        "severe_pairs_with_exact_zero_fraction": float(
            (severe & exact_zero).sum() / severe.sum()
        ),
        "published_tpm_ge_10_reanalysis_zero_count": int(
            published_high_reanalysis_zero.sum()
        ),
        "reanalysis_tpm_ge_10_published_zero_count": int(
            reanalysis_high_published_zero.sum()
        ),
        "genes_with_abs_log2_error_gt_2_in_any_sample": int(
            (severe_samples_per_gene > 0).sum()
        ),
        "genes_with_abs_log2_error_gt_2_in_all_samples": int(
            (severe_samples_per_gene == severe.shape[1]).sum()
        ),
        "genes_with_abs_log2_error_gt_2_in_at_least_30_samples": int(
            (severe_samples_per_gene >= 30).sum()
        ),
    }


def zero_transition_analysis(
    published_tpm: np.ndarray,
    reanalysis_tpm: np.ndarray,
    samples: list[str],
    gene_map: pd.DataFrame,
) -> tuple[go.Figure, dict[str, object], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    published_zero = published_tpm == 0
    reanalysis_zero = reanalysis_tpm == 0
    gained_any = published_zero & (reanalysis_tpm > 0)
    lost_any = (published_tpm > 0) & reanalysis_zero
    thresholds = [0.0, 0.1, 1.0, 5.0, 10.0]
    threshold_rows = []
    for threshold in thresholds:
        gained = published_zero & (
            (reanalysis_tpm > 0) if threshold == 0 else (reanalysis_tpm >= threshold)
        )
        lost = reanalysis_zero & (
            (published_tpm > 0) if threshold == 0 else (published_tpm >= threshold)
        )
        label = ">0" if threshold == 0 else f"≥{threshold:g}"
        threshold_rows.extend(
            [
                {
                    "nonzero_tpm_threshold": threshold,
                    "threshold_label": label,
                    "direction": "published 0 → reanalysis nonzero",
                    "gene_sample_pairs": int(gained.sum()),
                },
                {
                    "nonzero_tpm_threshold": threshold,
                    "threshold_label": label,
                    "direction": "published nonzero → reanalysis 0",
                    "gene_sample_pairs": int(lost.sum()),
                },
            ]
        )
    threshold_table = pd.DataFrame(threshold_rows)

    per_sample_rows = []
    for index, sample in enumerate(samples):
        per_sample_rows.append(
            {
                "sample": sample,
                "published_zero_to_reanalysis_nonzero": int(gained_any[:, index].sum()),
                "published_nonzero_to_reanalysis_zero": int(lost_any[:, index].sum()),
                "published_zero_to_reanalysis_ge_1": int(
                    (published_zero[:, index] & (reanalysis_tpm[:, index] >= 1)).sum()
                ),
                "published_ge_1_to_reanalysis_zero": int(
                    ((published_tpm[:, index] >= 1) & reanalysis_zero[:, index]).sum()
                ),
                "published_zero_to_reanalysis_ge_10": int(
                    (published_zero[:, index] & (reanalysis_tpm[:, index] >= 10)).sum()
                ),
                "published_ge_10_to_reanalysis_zero": int(
                    ((published_tpm[:, index] >= 10) & reanalysis_zero[:, index]).sum()
                ),
            }
        )
    per_sample = pd.DataFrame(per_sample_rows)

    gained_count = gained_any.sum(axis=1)
    lost_count = lost_any.sum(axis=1)
    gained_sum = np.where(gained_any, reanalysis_tpm, 0).sum(axis=1)
    lost_sum = np.where(lost_any, published_tpm, 0).sum(axis=1)
    gained_max = np.where(gained_any, reanalysis_tpm, -np.inf).max(axis=1)
    lost_max = np.where(lost_any, published_tpm, -np.inf).max(axis=1)
    gained_max[~np.isfinite(gained_max)] = np.nan
    lost_max[~np.isfinite(lost_max)] = np.nan
    per_gene = gene_map[
        ["published_id", "published_symbol", "reanalysis_id", "reanalysis_symbol"]
    ].copy()
    per_gene["published_zero_to_reanalysis_nonzero_samples"] = gained_count
    per_gene["published_zero_to_reanalysis_nonzero_mean_tpm"] = np.divide(
        gained_sum,
        gained_count,
        out=np.full(len(gained_count), np.nan),
        where=gained_count > 0,
    )
    per_gene["published_zero_to_reanalysis_nonzero_max_tpm"] = gained_max
    per_gene["published_nonzero_to_reanalysis_zero_samples"] = lost_count
    per_gene["published_nonzero_to_reanalysis_zero_mean_tpm"] = np.divide(
        lost_sum,
        lost_count,
        out=np.full(len(lost_count), np.nan),
        where=lost_count > 0,
    )
    per_gene["published_nonzero_to_reanalysis_zero_max_tpm"] = lost_max

    figure = make_subplots(
        rows=2,
        cols=2,
        specs=[[{}, {}], [{"colspan": 2}, None]],
        subplot_titles=(
            "Transition counts versus nonzero-side TPM threshold",
            "Expression on the nonzero side of exact-zero transitions",
            "Transitions by sample (nonzero side TPM ≥1)",
        ),
        horizontal_spacing=0.12,
        vertical_spacing=0.18,
    )
    direction_colors = {
        "published 0 → reanalysis nonzero": "#059669",
        "published nonzero → reanalysis 0": "#dc2626",
    }
    for direction, group in threshold_table.groupby("direction", sort=False):
        figure.add_trace(
            go.Bar(
                x=group["threshold_label"],
                y=group["gene_sample_pairs"],
                name=direction,
                marker_color=direction_colors[direction],
                showlegend=False,
                hovertemplate=(
                    direction
                    + "<br>nonzero TPM %{x}<br>%{y:,} gene×sample pairs<extra></extra>"
                ),
            ),
            row=1,
            col=1,
        )

    gained_values = np.log2(reanalysis_tpm[gained_any] + 1)
    lost_values = np.log2(published_tpm[lost_any] + 1)
    distribution_max = float(
        np.quantile(np.concatenate([gained_values, lost_values]), 0.999)
    )
    edges = np.linspace(0, distribution_max, 100)
    centers = (edges[:-1] + edges[1:]) / 2
    for values, name, color in (
        (gained_values, "published 0 → reanalysis nonzero", "#059669"),
        (lost_values, "published nonzero → reanalysis 0", "#dc2626"),
    ):
        counts, _ = np.histogram(values, bins=edges)
        figure.add_trace(
            go.Scatter(
                x=centers,
                y=counts / counts.sum(),
                mode="lines",
                name=name,
                legendgroup=name,
                showlegend=False,
                line={"color": color, "width": 2},
                hovertemplate="log2(TPM+1) %{x:.2f}<br>fraction %{y:.3%}<extra></extra>",
            ),
            row=1,
            col=2,
        )

    figure.add_trace(
        go.Bar(
            x=per_sample["sample"],
            y=per_sample["published_zero_to_reanalysis_ge_1"],
            name="published 0 → reanalysis ≥1",
            marker_color="#059669",
            showlegend=False,
        ),
        row=2,
        col=1,
    )
    figure.add_trace(
        go.Bar(
            x=per_sample["sample"],
            y=per_sample["published_ge_1_to_reanalysis_zero"],
            name="published ≥1 → reanalysis 0",
            marker_color="#dc2626",
            showlegend=False,
        ),
        row=2,
        col=1,
    )

    figure.update_layout(
        height=980,
        title="Exact zero ↔ nonzero transitions",
        template="plotly_white",
        barmode="group",
        showlegend=False,
        margin={"l": 70, "r": 40, "t": 125, "b": 145},
    )
    figure.update_xaxes(title_text="TPM on the nonzero side", row=1, col=1)
    figure.update_yaxes(title_text="Gene×sample pairs", row=1, col=1)
    figure.update_xaxes(title_text="log2(TPM + 1)", row=1, col=2)
    figure.update_yaxes(title_text="Fraction of transitions", row=1, col=2)
    figure.update_xaxes(tickangle=55, row=2, col=1)
    figure.update_yaxes(title_text="Genes per sample", row=2, col=1)

    threshold_lookup = threshold_table.pivot(
        index="nonzero_tpm_threshold", columns="direction", values="gene_sample_pairs"
    )
    summary: dict[str, object] = {
        "published_zero_to_reanalysis_nonzero_count": int(gained_any.sum()),
        "published_nonzero_to_reanalysis_zero_count": int(lost_any.sum()),
        "published_zero_to_reanalysis_ge_1_count": int(
            threshold_lookup.at[1.0, "published 0 → reanalysis nonzero"]
        ),
        "published_ge_1_to_reanalysis_zero_count": int(
            threshold_lookup.at[1.0, "published nonzero → reanalysis 0"]
        ),
        "published_zero_to_reanalysis_ge_10_count": int(
            threshold_lookup.at[10.0, "published 0 → reanalysis nonzero"]
        ),
        "published_ge_10_to_reanalysis_zero_count": int(
            threshold_lookup.at[10.0, "published nonzero → reanalysis 0"]
        ),
        "genes_ever_published_zero_to_reanalysis_nonzero": int((gained_count > 0).sum()),
        "genes_ever_published_nonzero_to_reanalysis_zero": int((lost_count > 0).sum()),
    }
    return figure, summary, threshold_table, per_sample, per_gene


def align_pca_axes(
    reference: np.ndarray,
    moving: np.ndarray,
    moving_variance: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    correlations = np.corrcoef(reference.T, moving.T)[:2, 2:]
    identity_score = abs(correlations[0, 0]) + abs(correlations[1, 1])
    swapped_score = abs(correlations[0, 1]) + abs(correlations[1, 0])
    order = np.array([0, 1]) if identity_score >= swapped_score else np.array([1, 0])
    aligned = moving[:, order].copy()
    variance = moving_variance[order].copy()
    for component in range(2):
        if np.corrcoef(reference[:, component], aligned[:, component])[0, 1] < 0:
            aligned[:, component] *= -1
    return aligned, variance


def select_variable_genes(
    transformed: np.ndarray,
    candidates: np.ndarray,
    top_genes: int,
) -> tuple[np.ndarray, str]:
    """Select the highest-variance genes without scaling gene variance."""
    variances = transformed.var(axis=1, ddof=1)
    return select_genes_by_variance(variances, candidates, top_genes)


def select_genes_by_variance(
    variances: np.ndarray,
    candidates: np.ndarray,
    top_genes: int,
) -> tuple[np.ndarray, str]:
    """Select candidate genes from a precomputed per-gene variance score."""
    candidate_variances = variances[candidates]
    variable = candidates[
        np.isfinite(candidate_variances) & (candidate_variances > 0)
    ]
    if len(variable) < 2:
        raise ValueError("Fewer than two variable matched genes remain for PCA")
    if top_genes <= 0 or top_genes >= len(variable):
        return variable, "all_variable_genes"
    order = np.argsort(variances[variable], kind="stable")
    return variable[order[-top_genes:]], "top_variable_genes"


def calculate_pca(
    published_log: np.ndarray,
    reanalysis_log: np.ndarray,
    min_mean_tpm: float,
    top_genes: int,
    published_tpm: np.ndarray,
    reanalysis_tpm: np.ndarray,
) -> tuple[dict[str, np.ndarray], dict[str, float | int]]:
    mean_tpm = (published_tpm.mean(axis=1) + reanalysis_tpm.mean(axis=1)) / 2
    candidates = np.flatnonzero(mean_tpm >= min_mean_tpm)
    if len(candidates) < 2:
        raise ValueError("Fewer than two matched genes remain for PCA")

    # Both comparison pages use the same procedure: log2(TPM + 1), independent
    # top-variable selection for each source, and PCA without per-gene variance
    # scaling. The joint panel makes its own selection on the combined profiles.
    published_keep, published_selection = select_variable_genes(
        published_log, candidates, top_genes
    )
    reanalysis_keep, reanalysis_selection = select_variable_genes(
        reanalysis_log, candidates, top_genes
    )
    average_within_source_variance = (
        published_log.var(axis=1, ddof=1)
        + reanalysis_log.var(axis=1, ddof=1)
    ) / 2
    joint_keep, joint_selection = select_genes_by_variance(
        average_within_source_variance, candidates, top_genes
    )

    published_x = published_log[published_keep].T
    reanalysis_x = reanalysis_log[reanalysis_keep].T

    published_model = PCA(n_components=2, random_state=42)
    reanalysis_model = PCA(n_components=2, random_state=42)
    published_scores = published_model.fit_transform(published_x)
    reanalysis_scores = reanalysis_model.fit_transform(reanalysis_x)
    reanalysis_scores, reanalysis_variance = align_pca_axes(
        published_scores,
        reanalysis_scores,
        reanalysis_model.explained_variance_ratio_,
    )

    published_proc, reanalysis_proc, disparity = procrustes(
        published_scores, reanalysis_scores
    )
    distance_pearson = safe_pearson(
        pdist(published_x), pdist(reanalysis_x)
    )
    distance_spearman = safe_spearman(
        pdist(published_x), pdist(reanalysis_x)
    )

    joint_x = np.vstack(
        [published_log[joint_keep].T, reanalysis_log[joint_keep].T]
    )
    joint_model = PCA(n_components=2, random_state=42)
    joint_scores = joint_model.fit_transform(joint_x)

    arrays = {
        "keep": joint_keep,
        "published_keep": published_keep,
        "reanalysis_keep": reanalysis_keep,
        "joint_keep": joint_keep,
        "published_scores": published_scores,
        "reanalysis_scores": reanalysis_scores,
        "published_procrustes": published_proc,
        "reanalysis_procrustes": reanalysis_proc,
        "joint_scores": joint_scores,
        "published_variance": published_model.explained_variance_ratio_,
        "reanalysis_variance": reanalysis_variance,
        "joint_variance": joint_model.explained_variance_ratio_,
    }
    metrics: dict[str, float | int] = {
        "genes_used": int(len(joint_keep)),
        "published_genes_used": int(len(published_keep)),
        "reanalysis_genes_used": int(len(reanalysis_keep)),
        "joint_genes_used": int(len(joint_keep)),
        "genes_passing_expression_filter": int(len(candidates)),
        "gene_selection": "source_specific_top_variance",
        "published_gene_selection": published_selection,
        "reanalysis_gene_selection": reanalysis_selection,
        "joint_gene_selection": joint_selection,
        "joint_selection_basis": "average_within_source_variance",
        "min_mean_tpm": float(min_mean_tpm),
        "transformation": "log2_tpm_plus_1",
        "per_gene_standardization": False,
        "protocol": "top_500_variable_log2_tpm_unscaled",
        "published_pc1_variance_pct": float(
            published_model.explained_variance_ratio_[0] * 100
        ),
        "published_pc2_variance_pct": float(
            published_model.explained_variance_ratio_[1] * 100
        ),
        "reanalysis_pc1_variance_pct": float(reanalysis_variance[0] * 100),
        "reanalysis_pc2_variance_pct": float(reanalysis_variance[1] * 100),
        "joint_pc1_variance_pct": float(joint_model.explained_variance_ratio_[0] * 100),
        "joint_pc2_variance_pct": float(joint_model.explained_variance_ratio_[1] * 100),
        "procrustes_disparity": float(disparity),
        "pairwise_sample_distance_pearson": distance_pearson,
        "pairwise_sample_distance_spearman": distance_spearman,
    }
    return arrays, metrics


def density_heatmap(
    x: np.ndarray,
    y: np.ndarray,
    x_range: tuple[float, float],
    y_range: tuple[float, float],
    bins: int = 100,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    counts, x_edges, y_edges = np.histogram2d(
        x, y, bins=bins, range=[x_range, y_range]
    )
    x_centers = (x_edges[:-1] + x_edges[1:]) / 2
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2
    return x_centers, y_centers, np.log10(counts.T + 1)


def error_figure(
    published_log: np.ndarray,
    reanalysis_log: np.ndarray,
    per_sample: pd.DataFrame,
) -> go.Figure:
    residual = (reanalysis_log - published_log).ravel()
    published_flat = published_log.ravel()
    reanalysis_flat = reanalysis_log.ravel()
    average = ((published_log + reanalysis_log) / 2).ravel()
    residual_range = tuple(np.quantile(residual, [0.001, 0.999]))
    expression_max = float(np.quantile(np.concatenate([published_flat, reanalysis_flat]), 0.999))

    figure = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=(
            "Published versus reanalysed expression",
            "Error distribution",
            "Error by sample (median and 10th–90th percentile)",
            "Error versus average expression",
        ),
        horizontal_spacing=0.12,
        vertical_spacing=0.16,
    )
    x, y, z = density_heatmap(
        published_flat,
        reanalysis_flat,
        (0.0, expression_max),
        (0.0, expression_max),
    )
    figure.add_trace(
        go.Heatmap(x=x, y=y, z=z, colorscale="Viridis", colorbar_title="log10(n+1)"),
        row=1,
        col=1,
    )
    figure.add_trace(
        go.Scatter(
            x=[0, expression_max],
            y=[0, expression_max],
            mode="lines",
            line={"color": "white", "dash": "dash", "width": 1},
            hoverinfo="skip",
            showlegend=False,
        ),
        row=1,
        col=1,
    )

    hist, edges = np.histogram(residual, bins=140, range=residual_range)
    centers = (edges[:-1] + edges[1:]) / 2
    figure.add_trace(
        go.Bar(x=centers, y=hist / hist.sum(), marker_color="#475569", showlegend=False),
        row=1,
        col=2,
    )

    medians = per_sample["median_signed_log2_error"].to_numpy()
    figure.add_trace(
        go.Scatter(
            x=per_sample["sample"],
            y=medians,
            mode="markers",
            marker={"color": "#7c3aed", "size": 7},
            error_y={
                "type": "data",
                "symmetric": False,
                "array": per_sample["error_q90"].to_numpy() - medians,
                "arrayminus": medians - per_sample["error_q10"].to_numpy(),
                "thickness": 1,
            },
            showlegend=False,
        ),
        row=2,
        col=1,
    )

    x, y, z = density_heatmap(
        average,
        residual,
        (0.0, expression_max),
        residual_range,
    )
    figure.add_trace(
        go.Heatmap(x=x, y=y, z=z, colorscale="Magma", showscale=False),
        row=2,
        col=2,
    )
    diagonal_x = np.linspace(0, expression_max, 300)
    for diagonal_y, label in (
        (2 * diagonal_x, "Published TPM = 0"),
        (-2 * diagonal_x, "Reanalysis TPM = 0"),
    ):
        visible = (diagonal_y >= residual_range[0]) & (diagonal_y <= residual_range[1])
        figure.add_trace(
            go.Scatter(
                x=diagonal_x[visible],
                y=diagonal_y[visible],
                mode="lines",
                line={"color": "#67e8f9", "dash": "dot", "width": 1},
                name=label,
                hovertemplate=label + "<extra></extra>",
                showlegend=False,
            ),
            row=2,
            col=2,
        )
    for row, col in ((1, 2), (2, 1), (2, 2)):
        figure.add_hline(y=0, line_dash="dash", line_color="#94a3b8", row=row, col=col)
    figure.update_xaxes(title_text="Published log2(TPM + 1)", row=1, col=1)
    figure.update_yaxes(title_text="Reanalysis log2(TPM + 1)", row=1, col=1)
    figure.update_xaxes(title_text="Reanalysis − published log2(TPM + 1)", row=1, col=2)
    figure.update_yaxes(title_text="Fraction of matched gene-sample pairs", row=1, col=2)
    figure.update_xaxes(tickangle=55, row=2, col=1)
    figure.update_yaxes(title_text="Reanalysis − published", row=2, col=1)
    figure.update_xaxes(title_text="Mean log2(TPM + 1)", row=2, col=2)
    figure.update_yaxes(title_text="Reanalysis − published", row=2, col=2)
    figure.update_layout(
        height=950,
        title="TPM agreement across matched genes and samples",
        template="plotly_white",
        margin={"l": 70, "r": 40, "t": 100, "b": 130},
    )
    return figure


def condition_palette(metadata: pd.DataFrame) -> dict[str, str]:
    conditions = list(dict.fromkeys(metadata["reproductive_state"]))
    return {
        condition: CONDITION_COLORS[index % len(CONDITION_COLORS)]
        for index, condition in enumerate(conditions)
    }


def add_condition_points(
    figure: go.Figure,
    scores: np.ndarray,
    metadata: pd.DataFrame,
    row: int,
    col: int,
    palette: dict[str, str],
    source: str,
    show_legend: bool,
    symbol: str = "circle",
) -> None:
    for condition, group in metadata.groupby("reproductive_state", sort=False):
        indices = group.index.to_numpy()
        figure.add_trace(
            go.Scatter(
                x=scores[indices, 0],
                y=scores[indices, 1],
                mode="markers",
                name=condition,
                legendgroup=condition,
                showlegend=show_legend,
                marker={
                    "color": palette[condition],
                    "symbol": symbol,
                    "size": 10,
                    "line": {"color": "white", "width": 0.8},
                },
                customdata=np.column_stack(
                    [group["sample"].to_numpy(), np.repeat(source, len(group))]
                ),
                hovertemplate=(
                    "%{customdata[0]}<br>%{customdata[1]}<br>"
                    + escape(condition)
                    + "<extra></extra>"
                ),
            ),
            row=row,
            col=col,
        )


def add_pair_lines(
    figure: go.Figure,
    first: np.ndarray,
    second: np.ndarray,
    row: int,
    col: int,
) -> None:
    x: list[float | None] = []
    y: list[float | None] = []
    for a, b in zip(first, second):
        x.extend([float(a[0]), float(b[0]), None])
        y.extend([float(a[1]), float(b[1]), None])
    figure.add_trace(
        go.Scatter(
            x=x,
            y=y,
            mode="lines",
            line={"color": "#cbd5e1", "width": 1},
            hoverinfo="skip",
            showlegend=False,
        ),
        row=row,
        col=col,
    )


def pca_figure(
    arrays: dict[str, np.ndarray],
    metadata: pd.DataFrame,
    group_label: str = "Reproductive state",
) -> go.Figure:
    pv = arrays["published_variance"] * 100
    rv = arrays["reanalysis_variance"] * 100
    jv = arrays["joint_variance"] * 100
    # Stacked, not side by side: three squeezed panels in one row made the
    # condition clusters unreadable.
    figure = make_subplots(
        rows=3,
        cols=1,
        subplot_titles=(
            f"Published PCA (PC1+PC2 {pv.sum():.1f}%)",
            f"Reanalysis PCA (PC1+PC2 {rv.sum():.1f}%)",
            f"Joint PCA of all {2 * len(metadata)} profiles (PC1+PC2 {jv.sum():.1f}%)",
        ),
        vertical_spacing=0.08,
    )
    palette = condition_palette(metadata)
    add_condition_points(
        figure,
        arrays["published_scores"],
        metadata,
        1,
        1,
        palette,
        "Published",
        True,
    )
    add_condition_points(
        figure,
        arrays["reanalysis_scores"],
        metadata,
        2,
        1,
        palette,
        "Reanalysis",
        False,
        "x",
    )

    sample_count = len(metadata)
    joint_published = arrays["joint_scores"][:sample_count]
    joint_reanalysis = arrays["joint_scores"][sample_count:]
    add_pair_lines(figure, joint_published, joint_reanalysis, 3, 1)
    add_condition_points(
        figure,
        joint_published,
        metadata,
        3,
        1,
        palette,
        "Published",
        False,
    )
    add_condition_points(
        figure,
        joint_reanalysis,
        metadata,
        3,
        1,
        palette,
        "Reanalysis",
        False,
        "x",
    )

    figure.update_xaxes(title_text=f"PC1 ({pv[0]:.1f}%)", row=1, col=1)
    figure.update_yaxes(title_text=f"PC2 ({pv[1]:.1f}%)", row=1, col=1)
    figure.update_xaxes(title_text=f"PC1 ({rv[0]:.1f}%)", row=2, col=1)
    figure.update_yaxes(title_text=f"PC2 ({rv[1]:.1f}%)", row=2, col=1)
    figure.update_xaxes(title_text=f"Joint PC1 ({jv[0]:.1f}%)", row=3, col=1)
    figure.update_yaxes(title_text=f"Joint PC2 ({jv[1]:.1f}%)", row=3, col=1)
    figure.update_layout(
        height=1750,
        title="Sample-level PCA comparison (joint PCA: circle = published, × = reanalysis)",
        template="plotly_white",
        legend={"title": group_label, "orientation": "v"},
        margin={"l": 65, "r": 285, "t": 105, "b": 65},
    )
    return figure


def sample_correlation_figure(
    published_log: np.ndarray,
    reanalysis_log: np.ndarray,
    samples: list[str],
) -> tuple[go.Figure, dict[str, float | int]]:
    matrix = np.corrcoef(published_log.T, reanalysis_log.T)[: len(samples), len(samples) :]
    best = np.argmax(matrix, axis=1)
    diagonal = np.diag(matrix)
    metrics: dict[str, float | int] = {
        "matching_sample_top1_count": int(np.sum(best == np.arange(len(samples)))),
        "matching_sample_top1_fraction": float(np.mean(best == np.arange(len(samples)))),
        "median_matching_sample_correlation": float(np.median(diagonal)),
        "min_matching_sample_correlation": float(np.min(diagonal)),
        "max_matching_sample_correlation": float(np.max(diagonal)),
    }
    figure = go.Figure(
        go.Heatmap(
            z=matrix,
            x=samples,
            y=samples,
            colorscale="Viridis",
            zmin=float(np.quantile(matrix, 0.01)),
            zmax=1.0,
            colorbar_title="Pearson r",
            hovertemplate="Published: %{y}<br>Reanalysis: %{x}<br>r=%{z:.4f}<extra></extra>",
        )
    )
    figure.update_layout(
        height=850,
        title="Cross-source sample correlation on matched log2(TPM + 1) genes",
        xaxis_title="Reanalysis sample",
        yaxis_title="Published sample",
        xaxis={"tickangle": 55},
        template="plotly_white",
        margin={"l": 170, "r": 70, "t": 90, "b": 180},
    )
    return figure, metrics


def pca_coordinates(
    arrays: dict[str, np.ndarray], metadata: pd.DataFrame
) -> pd.DataFrame:
    rows = []
    sample_count = len(metadata)
    embeddings = {
        "separate_pca": {
            "Published": arrays["published_scores"],
            "Reanalysis": arrays["reanalysis_scores"],
        },
        "procrustes_aligned_separate_pca": {
            "Published": arrays["published_procrustes"],
            "Reanalysis": arrays["reanalysis_procrustes"],
        },
        "joint_pca": {
            "Published": arrays["joint_scores"][:sample_count],
            "Reanalysis": arrays["joint_scores"][sample_count:],
        },
    }
    for embedding, by_source in embeddings.items():
        for source, scores in by_source.items():
            for index, row in metadata.iterrows():
                rows.append(
                    {
                        "embedding": embedding,
                        "source": source,
                        "sample": row["sample"],
                        "reproductive_state": row["reproductive_state"],
                        "component_1": scores[index, 0],
                        "component_2": scores[index, 1],
                    }
                )
    return pd.DataFrame(rows)


def zero_transition_tables_html(per_gene: pd.DataFrame, limit: int = 12) -> str:
    configurations = [
        (
            "Largest published 0 → reanalysis nonzero transitions",
            "published_zero_to_reanalysis_nonzero_samples",
            "published_zero_to_reanalysis_nonzero_mean_tpm",
            "published_zero_to_reanalysis_nonzero_max_tpm",
            "Reanalysis TPM",
        ),
        (
            "Largest published nonzero → reanalysis 0 transitions",
            "published_nonzero_to_reanalysis_zero_samples",
            "published_nonzero_to_reanalysis_zero_mean_tpm",
            "published_nonzero_to_reanalysis_zero_max_tpm",
            "Published TPM",
        ),
    ]
    tables = []
    for title, count_column, mean_column, max_column, tpm_label in configurations:
        selected = (
            per_gene.loc[per_gene[count_column] > 0]
            .sort_values([max_column, count_column], ascending=False, kind="stable")
            .head(limit)
        )
        body = "".join(
            "<tr>"
            f"<td>{escape(str(row.published_id))}</td>"
            f"<td>{escape(str(row.reanalysis_id))}</td>"
            f"<td>{int(getattr(row, count_column))}</td>"
            f"<td>{getattr(row, mean_column):,.2f}</td>"
            f"<td>{getattr(row, max_column):,.2f}</td>"
            "</tr>"
            for row in selected.itertuples(index=False)
        )
        tables.append(
            f"<div><h3>{escape(title)}</h3><table><thead><tr>"
            "<th>Published ID</th><th>Reanalysis ID</th><th>Samples</th>"
            f"<th>Mean {escape(tpm_label)}</th><th>Maximum {escape(tpm_label)}</th>"
            f"</tr></thead><tbody>{body}</tbody></table></div>"
        )
    return '<div class="table-grid">' + "".join(tables) + "</div>"


def render_report(
    output: Path,
    summary: dict[str, object],
    error_plot: go.Figure,
    zero_transition_plot: go.Figure,
    zero_transition_genes: pd.DataFrame,
    pca_plot: go.Figure,
    correlation_plot: go.Figure,
    report_title: str = "Published versus reanalysed ovary TPM",
    report_subtitle: str = (
        "Venkataraman et al. eLife 2023, 33 matched biological samples, "
        "log expression is <code>log2(TPM + 1)</code>."
    ),
) -> None:
    agreement = summary["agreement"]
    discordance = summary["discordance"]
    zero_transitions = summary["zero_transitions"]
    pca = summary["pca"]
    sample_identity = summary["sample_identity"]
    metrics = [
        ("Matched genes", f"{summary['matched_genes']:,} / {summary['published_genes']:,}"),
        ("Overall log-expression Pearson r", f"{agreement['pearson_log2_tpm_plus_1']:.4f}"),
        ("Median |log2(TPM+1) error|", f"{agreement['median_absolute_log2_error']:.3f}"),
        ("Within ±1 log2(TPM+1)", f"{agreement['fraction_abs_error_le_1']:.1%}"),
        (
            "Correct sample is top correlation",
            f"{sample_identity['matching_sample_top1_count']} / {summary['samples']}",
        ),
        ("Pairwise sample-distance Pearson r", f"{pca['pairwise_sample_distance_pearson']:.4f}"),
        ("PCA genes", f"{pca['genes_used']:,}"),
        (
            "Pairs with >2 log2 error",
            f"{discordance['abs_log2_error_gt_2_fraction']:.2%}",
        ),
        (
            "Severe pairs involving an exact zero",
            f"{discordance['severe_pairs_with_exact_zero_fraction']:.1%}",
        ),
    ]
    cards = "".join(
        f'<div class="card"><div class="label">{escape(label)}</div>'
        f'<div class="value">{escape(value)}</div></div>'
        for label, value in metrics
    )
    error_html = pio.to_html(error_plot, include_plotlyjs="inline", full_html=False)
    zero_transition_html = pio.to_html(
        zero_transition_plot, include_plotlyjs=False, full_html=False
    )
    zero_tables = zero_transition_tables_html(zero_transition_genes)
    pca_html = pio.to_html(pca_plot, include_plotlyjs=False, full_html=False)
    correlation_html = pio.to_html(
        correlation_plot, include_plotlyjs=False, full_html=False
    )
    pca_description = (
        "Each panel starts from one-to-one matched gene TPM transformed with "
        "<code>log2(TPM + 1)</code>. The separate panels retain their own "
        f"{pca['genes_used']:,} highest-variance genes. The joint panel ranks each "
        "gene by the average of its variance within the published and reanalysis "
        "profiles, so a processing-wide mean offset cannot by itself select a gene. "
        "PCA centers the retained gene columns but does not scale them to unit "
        "variance. This is not the paper's raw-count VST PCA."
    )
    output.write_text(
        f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{escape(report_title)}</title>
<style>
body {{ margin: 0; font-family: Inter, ui-sans-serif, system-ui, sans-serif; color: #0f172a; background: #f8fafc; }}
main {{ max-width: 1450px; margin: 0 auto; padding: 36px 28px 70px; }}
h1 {{ margin: 0 0 8px; font-size: 30px; }}
.subtitle {{ color: #475569; margin-bottom: 26px; line-height: 1.5; }}
.cards {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(210px, 1fr)); gap: 12px; margin-bottom: 24px; }}
.card, section {{ background: white; border: 1px solid #e2e8f0; border-radius: 12px; box-shadow: 0 1px 2px rgba(15,23,42,.04); }}
.card {{ padding: 16px; }} .label {{ color: #64748b; font-size: 13px; }} .value {{ font-size: 23px; font-weight: 650; margin-top: 5px; }}
section {{ padding: 14px; margin: 18px 0; overflow-x: auto; }}
.note {{ color: #475569; font-size: 14px; line-height: 1.55; margin: 10px 6px 2px; }}
.table-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(520px, 1fr)); gap: 22px; padding: 4px 12px 18px; }}
h3 {{ margin: 8px 0 10px; font-size: 16px; }}
table {{ width: 100%; border-collapse: collapse; font-size: 12px; }}
th, td {{ border-bottom: 1px solid #e2e8f0; padding: 7px 8px; text-align: right; }}
th:first-child, td:first-child, th:nth-child(2), td:nth-child(2) {{ text-align: left; }}
th {{ color: #475569; background: #f8fafc; position: sticky; top: 0; }}
code {{ background: #f1f5f9; padding: 2px 5px; border-radius: 5px; }}
</style></head><body><main>
<h1>{escape(report_title)}</h1>
<div class="subtitle">{report_subtitle}</div>
<div class="cards">{cards}</div>
<section>{error_html}<p class="note">Errors are reanalysis minus published values. The diagonal bands are forced by the coordinates: when published TPM is zero, error = +2 × average log-expression; when reanalysis TPM is zero, error = −2 × average log-expression. Of the {discordance['abs_log2_error_gt_2_count']:,} pairs with absolute error &gt;2, {discordance['severe_pairs_with_exact_zero_fraction']:.1%} contain an exact zero. Published TPM was ≥10 while reanalysis TPM was zero in {discordance['published_tpm_ge_10_reanalysis_zero_count']:,} pairs; the reverse occurred in {discordance['reanalysis_tpm_ge_10_published_zero_count']:,}. Density plots clip only the outer 0.1% for readable axes; summary metrics use the full distribution.</p></section>
<section>{zero_transition_html}<p class="note"><strong>Green</strong> denotes published 0 → reanalysis nonzero; <strong>red</strong> denotes published nonzero → reanalysis 0. An exact zero can mean no compatible fragments were assigned under that quantification model; it is not a universal biological absence threshold. Threshold bars therefore ask how large the value is on the nonzero side. There are {zero_transitions['published_zero_to_reanalysis_nonzero_count']:,} exact published 0 → reanalysis nonzero pairs and {zero_transitions['published_nonzero_to_reanalysis_zero_count']:,} published nonzero → reanalysis 0 pairs. At a nonzero-side threshold of 1 TPM these fall to {zero_transitions['published_zero_to_reanalysis_ge_1_count']:,} and {zero_transitions['published_ge_1_to_reanalysis_zero_count']:,}; at 10 TPM, {zero_transitions['published_zero_to_reanalysis_ge_10_count']:,} and {zero_transitions['published_ge_10_to_reanalysis_zero_count']:,}.</p>{zero_tables}</section>
<section>{pca_html}<p class="note">{pca_description} Separate PCAs compare biological geometry; the joint PCA also exposes method-specific shifts.</p></section>
<section>{correlation_html}<p class="note">The diagonal compares the same biological sample across processing methods. A diagonal maximum in each row argues against sample swaps.</p></section>
</main></body></html>""",
        encoding="utf-8",
    )


def render_zero_transition_report(
    output: Path,
    summary: dict[str, object],
    zero_transition_plot: go.Figure,
) -> None:
    zero_transitions = summary["zero_transitions"]
    plot_html = pio.to_html(
        zero_transition_plot, include_plotlyjs="inline", full_html=False
    )
    output.write_text(
        f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>eLife ovary zero-to-nonzero TPM transitions</title>
<style>
body {{ margin: 0; font-family: Inter, ui-sans-serif, system-ui, sans-serif; color: #0f172a; background: #f8fafc; }}
main {{ max-width: 1450px; margin: 0 auto; padding: 34px 28px 60px; }}
h1 {{ margin: 0 0 8px; font-size: 30px; }}
.subtitle {{ color: #475569; margin-bottom: 24px; line-height: 1.5; }}
section {{ background: white; border: 1px solid #e2e8f0; border-radius: 12px; box-shadow: 0 1px 2px rgba(15,23,42,.04); padding: 14px; overflow-x: auto; }}
.note {{ color: #475569; font-size: 14px; line-height: 1.6; margin: 8px 12px 14px; }}
</style></head><body><main>
<h1>Exact zero ↔ nonzero TPM transitions</h1>
<div class="subtitle">Published Venkataraman et al. eLife 2023 TPM versus the independent Salmon reanalysis of the same 33 ovary samples.</div>
<section>{plot_html}<p class="note"><strong>Green:</strong> published 0 → reanalysis nonzero. <strong>Red:</strong> published nonzero → reanalysis 0. Across all matched gene×sample pairs, the two directions contain {zero_transitions['published_zero_to_reanalysis_nonzero_count']:,} and {zero_transitions['published_nonzero_to_reanalysis_zero_count']:,} exact transitions. At a nonzero-side threshold of 1 TPM, they contain {zero_transitions['published_zero_to_reanalysis_ge_1_count']:,} and {zero_transitions['published_ge_1_to_reanalysis_zero_count']:,} transitions.</p></section>
</main></body></html>""",
        encoding="utf-8",
    )


def zero_transition_tables_records(
    per_gene: pd.DataFrame, limit: int = 12
) -> list[dict[str, object]]:
    """The same top-N transition tables the HTML report shows, as plain records."""
    configurations = [
        (
            "Largest published 0 to reanalysis nonzero transitions",
            "published_zero_to_reanalysis_nonzero_samples",
            "published_zero_to_reanalysis_nonzero_mean_tpm",
            "published_zero_to_reanalysis_nonzero_max_tpm",
            "Reanalysis TPM",
        ),
        (
            "Largest published nonzero to reanalysis 0 transitions",
            "published_nonzero_to_reanalysis_zero_samples",
            "published_nonzero_to_reanalysis_zero_mean_tpm",
            "published_nonzero_to_reanalysis_zero_max_tpm",
            "Published TPM",
        ),
    ]
    tables = []
    for title, count_column, mean_column, max_column, tpm_label in configurations:
        selected = (
            per_gene.loc[per_gene[count_column] > 0]
            .sort_values([max_column, count_column], ascending=False, kind="stable")
            .head(limit)
        )
        tables.append(
            {
                "title": title,
                "tpm_label": tpm_label,
                "rows": [
                    {
                        "Published ID": str(row.published_id),
                        "Reanalysis ID": str(row.reanalysis_id),
                        "Samples": int(getattr(row, count_column)),
                        f"Mean {tpm_label}": float(getattr(row, mean_column)),
                        f"Maximum {tpm_label}": float(getattr(row, max_column)),
                    }
                    for row in selected.itertuples(index=False)
                ],
            }
        )
    return tables


def export_figure_bundle(
    output: Path,
    summary: dict[str, object],
    figures: dict[str, go.Figure],
    zero_transition_genes: pd.DataFrame,
) -> None:
    """Write the figures as plotly JSON so the atlas can render them natively.

    The standalone HTML report stays for sharing, but embedding it in the app
    means an iframe with its own scrollbar and its own plotly copy. Emitting the
    figures lets the atlas draw them as ordinary page content with its own
    theme.
    """
    def plain(value: object) -> object:
        """Recursively cast numpy arrays/scalars to plain JSON types.

        pio.to_json encodes numpy arrays as base64 "bdata". Rebuilding a figure
        from that in the app mis-scaled the density heatmaps (a Magma panel that
        should be black rendered bright purple), so the bundle carries plain
        numbers instead.
        """
        if isinstance(value, np.ndarray):
            return [plain(item) for item in value.tolist()]
        if isinstance(value, dict):
            if {"dtype", "bdata"} <= set(value):
                decoded = np.frombuffer(
                    base64.b64decode(value["bdata"]), dtype=value["dtype"]
                )
                if value.get("shape"):
                    shape = tuple(int(n) for n in str(value["shape"]).split(","))
                    decoded = decoded.reshape(shape)
                return decoded.tolist()
            return {key: plain(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [plain(item) for item in value]
        if isinstance(value, (np.integer, np.floating)):
            return value.item()
        return value

    bundle = {
        "summary": summary,
        "figures": {
            name: plain(figure.to_plotly_json()) for name, figure in figures.items()
        },
        "zero_transition_tables": zero_transition_tables_records(zero_transition_genes),
    }
    output.write_text(json.dumps(bundle), encoding="utf-8")


def main() -> None:
    root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--published",
        type=Path,
        default=root / "expression/elife_80489_tpm.tsv.gz",
    )
    parser.add_argument(
        "--reanalysis",
        type=Path,
        default=root / "expression/elife_80489_salmon_gene_tpm.tsv.gz",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=root / "expression/elife_80489_samples.tsv",
    )
    parser.add_argument("--crosswalk", type=Path, required=True)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "analysis/results/elife_tpm_comparison",
    )
    parser.add_argument(
        "--top-variable-genes",
        type=int,
        default=500,
        help=(
            "Number of highest-variance genes selected independently for each PCA; "
            "0 uses all variable matched genes"
        ),
    )
    parser.add_argument("--min-mean-tpm", type=float, default=0.0)
    parser.add_argument("--group-label", default="Reproductive state")
    parser.add_argument(
        "--report-title", default="Published versus reanalysed ovary TPM"
    )
    parser.add_argument(
        "--report-subtitle",
        default=(
            "Venkataraman et al. eLife 2023, 33 matched biological samples, "
            "log expression is <code>log2(TPM + 1)</code>."
        ),
    )
    parser.add_argument(
        "--report-filename", default="elife_ovary_tpm_full_report.html"
    )
    args = parser.parse_args()

    published, published_samples = load_matrix(args.published)
    reanalysis, reanalysis_samples = load_matrix(args.reanalysis)
    if published_samples != reanalysis_samples:
        raise ValueError("Published and reanalysis sample columns are not identical")
    metadata = pd.read_csv(args.metadata, sep="\t", dtype=str).fillna("")
    metadata = metadata.set_index("sample").loc[published_samples].reset_index()
    crosswalk = pd.read_csv(args.crosswalk, sep="\t", dtype=str).fillna("")
    gene_map, mapping_diagnostics = build_gene_map(
        published, reanalysis, crosswalk
    )

    published_lookup = published.set_index("IDs")
    reanalysis_lookup = reanalysis.set_index("IDs")
    published_tpm = published_lookup.loc[
        gene_map["published_id"], published_samples
    ].to_numpy(dtype=float)
    reanalysis_tpm = reanalysis_lookup.loc[
        gene_map["reanalysis_id"], published_samples
    ].to_numpy(dtype=float)
    published_log = np.log2(published_tpm + 1.0)
    reanalysis_log = np.log2(reanalysis_tpm + 1.0)

    per_sample, agreement = error_metrics(
        published_log, reanalysis_log, published_samples
    )
    discordance = discordance_metrics(
        published_tpm, reanalysis_tpm, published_log, reanalysis_log
    )
    (
        zero_transition_plot,
        zero_transition_summary,
        zero_transition_thresholds,
        zero_transition_samples,
        zero_transition_genes,
    ) = zero_transition_analysis(
        published_tpm,
        reanalysis_tpm,
        published_samples,
        gene_map,
    )
    arrays, pca_metrics = calculate_pca(
        published_log,
        reanalysis_log,
        args.min_mean_tpm,
        args.top_variable_genes,
        published_tpm,
        reanalysis_tpm,
    )
    correlation_plot, sample_identity = sample_correlation_figure(
        published_log, reanalysis_log, published_samples
    )

    gene_residual = reanalysis_log - published_log
    gene_map["published_mean_tpm"] = published_tpm.mean(axis=1)
    gene_map["reanalysis_mean_tpm"] = reanalysis_tpm.mean(axis=1)
    gene_map["median_signed_log2_error"] = np.median(gene_residual, axis=1)
    gene_map["median_absolute_log2_error"] = np.median(
        np.abs(gene_residual), axis=1
    )
    gene_map["sample_profile_pearson"] = [
        safe_pearson(published_log[index], reanalysis_log[index])
        for index in range(len(gene_map))
    ]
    gene_map["used_for_published_pca"] = False
    gene_map["used_for_reanalysis_pca"] = False
    gene_map["used_for_joint_pca"] = False
    gene_map.loc[arrays["published_keep"], "used_for_published_pca"] = True
    gene_map.loc[arrays["reanalysis_keep"], "used_for_reanalysis_pca"] = True
    gene_map.loc[arrays["joint_keep"], "used_for_joint_pca"] = True
    gene_map["used_for_pca"] = gene_map[
        [
            "used_for_published_pca",
            "used_for_reanalysis_pca",
            "used_for_joint_pca",
        ]
    ].any(axis=1)

    summary: dict[str, object] = {
        "published_genes": int(len(published)),
        "reanalysis_genes": int(len(reanalysis)),
        "matched_genes": int(len(gene_map)),
        "matched_gene_fraction_of_published": float(len(gene_map) / len(published)),
        "samples": int(len(published_samples)),
        "mapping_methods": {
            key: int(value)
            for key, value in gene_map["mapping_method"].value_counts().items()
        },
        "mapping_diagnostics": mapping_diagnostics,
        "published_sample_tpm_sum_min": float(
            published[published_samples].sum().min()
        ),
        "published_sample_tpm_sum_max": float(
            published[published_samples].sum().max()
        ),
        "reanalysis_sample_tpm_sum_min": float(
            reanalysis[reanalysis_samples].sum().min()
        ),
        "reanalysis_sample_tpm_sum_max": float(
            reanalysis[reanalysis_samples].sum().max()
        ),
        "agreement": agreement,
        "discordance": discordance,
        "zero_transitions": zero_transition_summary,
        "sample_identity": sample_identity,
        "pca": pca_metrics,
    }

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    gene_map.to_csv(
        output_dir / "matched_genes.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    per_sample.merge(
        metadata[["sample", "reproductive_state", "replicate"]],
        on="sample",
        how="left",
        validate="one_to_one",
    ).to_csv(output_dir / "per_sample_metrics.tsv", sep="\t", index=False)
    pca_coordinates(arrays, metadata).to_csv(
        output_dir / "pca_coordinates.tsv", sep="\t", index=False
    )
    zero_transition_thresholds.to_csv(
        output_dir / "zero_transition_thresholds.tsv", sep="\t", index=False
    )
    zero_transition_samples.merge(
        metadata[["sample", "reproductive_state", "replicate"]],
        on="sample",
        how="left",
        validate="one_to_one",
    ).to_csv(output_dir / "zero_transition_per_sample.tsv", sep="\t", index=False)
    zero_transition_genes.to_csv(
        output_dir / "zero_transition_genes.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    error_plot = error_figure(published_log, reanalysis_log, per_sample)
    pca_plot = pca_figure(arrays, metadata, args.group_label)
    render_report(
        output_dir / "report.html",
        summary,
        error_plot,
        zero_transition_plot,
        zero_transition_genes,
        pca_plot,
        correlation_plot,
        args.report_title,
        args.report_subtitle,
    )
    export_figure_bundle(
        output_dir / "figures.json",
        summary,
        {
            "error": error_plot,
            "zero_transition": zero_transition_plot,
            "pca": pca_plot,
            "correlation": correlation_plot,
        },
        zero_transition_genes,
    )
    sendable = output_dir / args.report_filename
    sendable.write_bytes((output_dir / "report.html").read_bytes())
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"Wrote sendable report to {sendable}")


if __name__ == "__main__":
    main()
