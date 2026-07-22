"""Load, harmonize, search, and summarize the published TPM matrices."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from typing import Iterable

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class ExpressionDataset:
    key: str
    label: str
    paper: str
    annotation_version: str
    genes: pd.DataFrame
    values: pd.DataFrame
    samples: pd.DataFrame

    @property
    def sample_columns(self) -> list[str]:
        return list(self.values.columns)


DATASET_ORDER = ("elife", "neuro_ru", "neuro_legacy")

CONDITION_LABELS = {
    "BF": "Blood-fed",
    "O": "Gravid / oviposition-stage",
    "SF": "Non-blood-fed / sugar-fed",
    "male": "Male",
}

PAPER_FAMILY_LABELS = {
    "IR": "Ionotropic receptors (IR)",
    "OR": "Odorant receptors (OR)",
    "GR": "Gustatory receptors (GR)",
    "OBP": "Odorant-binding proteins (OBP)",
}

PAPER_ANNOTATION_COLUMNS = [
    "paper_gene_family",
    "orthodb_category",
    "drosophila_ortholog",
    "drosophila_blastx_hits",
    "naming_evidence",
]


def _clean(value: object) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    return str(value).strip()


def normalize_alias(value: object) -> str:
    """Normalize an identifier for forgiving exact matching."""
    text = _clean(value).casefold()
    return re.sub(r"[^a-z0-9]", "", text)


def canonical_symbol(symbol: object, stable_id: object = "") -> str:
    """Normalize known historical naming conventions without discarding IDs."""
    raw = _clean(symbol)
    stable = _clean(stable_id).upper()
    normalized = normalize_alias(raw)
    if stable == "AAEL005776" or normalized in {"orco", "aaegor7"}:
        return "Orco"
    if raw.casefold().startswith("aaeg") and len(raw) > 4:
        raw = raw[4:]
    return raw or _clean(stable_id)


def classify_family(symbol: object) -> str:
    normalized = normalize_alias(symbol).upper()
    if normalized == "ORCO" or re.fullmatch(r"OR\d+[A-Z]*", normalized):
        return "Odorant receptors (OR)"
    if re.fullmatch(r"IR\d+[A-Z]*", normalized):
        return "Ionotropic receptors (IR)"
    if re.fullmatch(r"GR\d+[A-Z]*", normalized):
        return "Gustatory receptors (GR)"
    if re.fullmatch(r"OBP\d+[A-Z]*", normalized):
        return "Odorant-binding proteins (OBP)"
    return "Other"


def _finalize_genes(
    annotations: pd.DataFrame,
    stable_id: pd.Series,
    internal_id: pd.Series,
    raw_symbol: pd.Series,
) -> pd.DataFrame:
    genes = annotations.copy().reset_index(drop=True)
    genes.insert(0, "row_id", np.arange(len(genes), dtype=int))
    genes["stable_id"] = stable_id.map(_clean).to_numpy()
    genes["internal_id"] = internal_id.map(_clean).to_numpy()
    genes["raw_symbol"] = raw_symbol.map(_clean).to_numpy()
    genes["canonical_symbol"] = [
        canonical_symbol(symbol, stable)
        for symbol, stable in zip(genes["raw_symbol"], genes["stable_id"])
    ]
    genes["display_name"] = genes["canonical_symbol"].where(
        genes["canonical_symbol"].ne(""), genes["stable_id"]
    )
    genes["family"] = genes["canonical_symbol"].map(classify_family)
    for column in PAPER_ANNOTATION_COLUMNS:
        genes[column] = ""

    def aliases(row: pd.Series) -> tuple[str, ...]:
        values = {
            row["stable_id"],
            row["internal_id"],
            row["raw_symbol"],
            row["canonical_symbol"],
        }
        if row["canonical_symbol"].casefold() == "orco":
            values.update({"Orco", "AaegOr7", "Or7", "AAEL005776"})
        return tuple(sorted(value for value in values if value))

    genes["aliases"] = genes.apply(aliases, axis=1)
    genes["search_text"] = genes["aliases"].map(
        lambda values: " | ".join(values)
    )
    genes["search_normalized"] = genes["aliases"].map(
        lambda values: " | ".join(normalize_alias(value) for value in values)
    )
    return genes


def _add_paper_annotations(genes: pd.DataFrame, path: Path) -> pd.DataFrame:
    annotations = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    annotations = annotations.drop_duplicates("vectorbase_identifier").set_index(
        "vectorbase_identifier"
    )
    enriched = genes.copy()
    for column in PAPER_ANNOTATION_COLUMNS:
        enriched[column] = enriched["stable_id"].map(annotations[column]).fillna("")
    paper_family = enriched["paper_gene_family"].map(PAPER_FAMILY_LABELS)
    enriched.loc[paper_family.notna(), "family"] = paper_family.dropna()
    return enriched


def _load_bmc_samples(path: Path, sample_columns: Iterable[str]) -> pd.DataFrame:
    samples = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    samples = samples[samples["library_id"].isin(sample_columns)].copy()
    samples["sample"] = samples["library_id"]
    samples["condition_label"] = samples["condition"].map(CONDITION_LABELS).fillna(
        samples["condition"]
    )
    samples["tissue"] = samples["tissue"].str.replace(
        "abdominaltip", "abdominal tip", regex=False
    )
    samples["tissue_condition"] = (
        samples["tissue"].str.title() + " · " + samples["condition_label"]
    )
    samples["reproductive_state"] = samples["condition_label"]
    samples = samples.set_index("sample", drop=False).reindex(list(sample_columns))
    return samples


def _load_elife_samples(path: Path, sample_columns: Iterable[str]) -> pd.DataFrame:
    samples = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    samples["sample"] = samples["sample"].astype(str)
    samples["condition"] = samples["reproductive_state"]
    samples["condition_label"] = samples["reproductive_state"]
    samples["tissue"] = "ovary"
    samples["sex"] = "female"
    samples["tissue_condition"] = "Ovary · " + samples["reproductive_state"]
    samples = samples.set_index("sample", drop=False).reindex(list(sample_columns))
    return samples


def _generic_samples(sample_columns: Iterable[str]) -> pd.DataFrame:
    samples = pd.DataFrame({"sample": list(sample_columns)})
    samples["condition"] = ""
    samples["condition_label"] = "Unspecified"
    samples["reproductive_state"] = "Unspecified"
    samples["tissue"] = "unspecified"
    samples["sex"] = ""
    samples["tissue_condition"] = "Unspecified"
    return samples.set_index("sample", drop=False)


def _read_matrix(path: Path, annotation_columns: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    frame = pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
    annotations = frame[annotation_columns].fillna("")
    values = frame.drop(columns=annotation_columns).apply(pd.to_numeric, errors="coerce").fillna(0.0)
    values.index = np.arange(len(values), dtype=int)
    return annotations, values


def load_nfcore_dataset(
    path: Path | str,
    key: str,
    symbol_crosswalk: dict[str, str] | None = None,
) -> ExpressionDataset:
    """Load an nf-core/rnaseq merged gene-TPM table.

    Supported examples include salmon.merged.gene_tpm.tsv and
    rsem.merged.gene_tpm.tsv (optionally gzip-compressed).
    """
    path = Path(path)
    frame = pd.read_csv(path, sep="\t", compression="infer", low_memory=False)
    if frame.shape[1] < 2:
        raise ValueError(f"nf-core TPM matrix needs a gene column and samples: {path}")
    gene_column = frame.columns[0]
    stable_ids = frame[gene_column].map(_clean)
    crosswalk = symbol_crosswalk or {}
    symbols = stable_ids.map(crosswalk).fillna("")
    annotations = pd.DataFrame({"gene_id": stable_ids})
    genes = _finalize_genes(
        annotations,
        stable_ids,
        pd.Series([""] * len(frame)),
        symbols,
    )
    values = frame.drop(columns=[gene_column]).apply(pd.to_numeric, errors="coerce")
    if values.isna().all(axis=None):
        raise ValueError(f"No numeric TPM values found in {path}")
    values = values.fillna(0.0)
    values.index = np.arange(len(values), dtype=int)
    label = path.name
    for suffix in (".gz", ".tsv"):
        if label.endswith(suffix):
            label = label[: -len(suffix)]
    return ExpressionDataset(
        key=key,
        label=f"nf-core · {label}",
        paper="Local nf-core/rnaseq import",
        annotation_version="Identifiers from imported matrix",
        genes=genes,
        values=values,
        samples=_generic_samples(values.columns),
    )


def load_datasets(expression_dir: Path | str) -> dict[str, ExpressionDataset]:
    expression_dir = Path(expression_dir)

    ru_annotations, ru_values = _read_matrix(
        expression_dir / "neurotranscriptome_2016_aaegl_ru_tpm.tsv.gz",
        ["Vectorbase Identifier", "Internal gene ID", "Display name"],
    )
    ru_genes = _finalize_genes(
        ru_annotations,
        ru_annotations["Vectorbase Identifier"],
        ru_annotations["Internal gene ID"],
        ru_annotations["Display name"],
    )
    annotation_path = expression_dir / "neurotranscriptome_2016_gene_annotations.tsv"
    ru_genes = _add_paper_annotations(ru_genes, annotation_path)
    crosswalk = {
        stable: symbol
        for stable, symbol in zip(ru_genes["stable_id"], ru_genes["raw_symbol"])
        if stable and symbol
    }
    bmc_samples = _load_bmc_samples(
        expression_dir / "neurotranscriptome_2016_samples.tsv", ru_values.columns
    )

    legacy_annotations, legacy_values = _read_matrix(
        expression_dir / "neurotranscriptome_2016_aaegl_3_3_tpm.tsv.gz",
        ["gene"],
    )
    legacy_symbols = legacy_annotations["gene"].map(crosswalk).fillna("")
    legacy_genes = _finalize_genes(
        legacy_annotations,
        legacy_annotations["gene"],
        pd.Series([""] * len(legacy_annotations)),
        legacy_symbols,
    )
    legacy_genes = _add_paper_annotations(legacy_genes, annotation_path)
    legacy_samples = _load_bmc_samples(
        expression_dir / "neurotranscriptome_2016_samples.tsv", legacy_values.columns
    )

    elife_annotations, elife_values = _read_matrix(
        expression_dir / "elife_80489_tpm.tsv.gz",
        ["IDs", "Symbols"],
    )
    elife_genes = _finalize_genes(
        elife_annotations,
        elife_annotations["IDs"],
        pd.Series([""] * len(elife_annotations)),
        elife_annotations["Symbols"],
    )
    elife_samples = _load_elife_samples(
        expression_dir / "elife_80489_samples.tsv", elife_values.columns
    )

    datasets = {
        "elife": ExpressionDataset(
            key="elife",
            label="Drought resilience · ovary time course",
            paper="Venkataraman et al., eLife 2023",
            annotation_version="Published eLife gene symbols",
            genes=elife_genes,
            values=elife_values,
            samples=elife_samples,
        ),
        "neuro_ru": ExpressionDataset(
            key="neuro_ru",
            label="Neurotranscriptome · updated AaegL.RU",
            paper="Matthews et al., BMC Genomics 2016",
            annotation_version="AaegL.RU (recommended)",
            genes=ru_genes,
            values=ru_values,
            samples=bmc_samples,
        ),
        "neuro_legacy": ExpressionDataset(
            key="neuro_legacy",
            label="Neurotranscriptome · legacy AaegL3.3",
            paper="Matthews et al., BMC Genomics 2016",
            annotation_version="AaegL3.3 (compatibility)",
            genes=legacy_genes,
            values=legacy_values,
            samples=legacy_samples,
        ),
    }

    import_dir = expression_dir / "imports"
    if import_dir.exists():
        imported_paths = sorted(import_dir.glob("*gene_tpm.tsv")) + sorted(
            import_dir.glob("*gene_tpm.tsv.gz")
        )
        used_keys: set[str] = set(datasets)
        for position, path in enumerate(imported_paths, start=1):
            slug = re.sub(r"[^a-z0-9]+", "_", path.name.casefold()).strip("_")
            key = f"nfcore_{slug or position}"
            while key in used_keys:
                key = f"{key}_{position}"
            used_keys.add(key)
            datasets[key] = load_nfcore_dataset(path, key, crosswalk)
    return datasets


def search_genes(dataset: ExpressionDataset, query: str, mode: str = "exact") -> pd.DataFrame:
    query = query.strip()
    if not query:
        return dataset.genes.iloc[0:0].copy()
    normalized = normalize_alias(query)
    if mode == "exact":
        mask = dataset.genes["aliases"].map(
            lambda values: normalized in {normalize_alias(value) for value in values}
        )
    elif mode == "prefix":
        mask = dataset.genes["aliases"].map(
            lambda values: any(normalize_alias(value).startswith(normalized) for value in values)
        )
    elif mode == "contains":
        mask = dataset.genes["search_normalized"].str.contains(
            re.escape(normalized), case=False, regex=True
        )
    elif mode == "regex":
        try:
            pattern = re.compile(query, re.IGNORECASE)
        except re.error as exc:
            raise ValueError(f"Invalid regular expression: {exc}") from exc
        mask = dataset.genes["search_text"].map(lambda value: bool(pattern.search(value)))
    else:
        raise ValueError(f"Unknown search mode: {mode}")
    return dataset.genes.loc[mask].copy()


def family_members(dataset: ExpressionDataset, family: str) -> pd.DataFrame:
    return dataset.genes.loc[dataset.genes["family"].eq(family)].copy()


def expression_long(dataset: ExpressionDataset, genes: pd.DataFrame) -> pd.DataFrame:
    if genes.empty:
        return pd.DataFrame()
    row_ids = genes["row_id"].astype(int).tolist()
    labels = genes.set_index("row_id")["display_name"].to_dict()
    stable_ids = genes.set_index("row_id")["stable_id"].to_dict()
    subset = dataset.values.loc[row_ids].copy()
    subset.insert(0, "row_id", row_ids)
    long = subset.melt(id_vars="row_id", var_name="sample", value_name="tpm")
    long["gene"] = long["row_id"].map(labels)
    long["stable_id"] = long["row_id"].map(stable_ids)
    metadata = dataset.samples.reset_index(drop=True)
    long = long.merge(metadata, on="sample", how="left", suffixes=("", "_sample"))
    long["dataset"] = dataset.label
    return long


def gene_statistics(dataset: ExpressionDataset, genes: pd.DataFrame) -> pd.DataFrame:
    if genes.empty:
        return pd.DataFrame()
    rows = []
    for gene in genes.itertuples(index=False):
        expression = dataset.values.loc[int(gene.row_id)]
        max_sample = str(expression.idxmax())
        sample_meta = dataset.samples.loc[max_sample]
        rows.append(
            {
                "gene": gene.display_name,
                "stable_id": gene.stable_id,
                "raw_symbol": gene.raw_symbol,
                "family": gene.family,
                "paper_gene_family": gene.paper_gene_family,
                "orthodb_category": gene.orthodb_category,
                "drosophila_ortholog": gene.drosophila_ortholog,
                "drosophila_blastx_hits": gene.drosophila_blastx_hits,
                "naming_evidence": gene.naming_evidence,
                "mean_tpm": float(expression.mean()),
                "median_tpm": float(expression.median()),
                "max_tpm": float(expression.max()),
                "detected_pct": float((expression >= 1.0).mean() * 100.0),
                "top_sample": max_sample,
                "top_context": _clean(sample_meta.get("tissue_condition", "")),
            }
        )
    return pd.DataFrame(rows)


def matrix_for_genes(dataset: ExpressionDataset, genes: pd.DataFrame) -> pd.DataFrame:
    """Return annotations plus the original per-sample TPM columns."""
    if genes.empty:
        return pd.DataFrame()
    row_ids = genes["row_id"].astype(int).tolist()
    annotations = genes.set_index("row_id").loc[row_ids, [
        "display_name",
        "stable_id",
        "internal_id",
        "raw_symbol",
        "family",
        "paper_gene_family",
        "orthodb_category",
        "drosophila_ortholog",
        "drosophila_blastx_hits",
        "naming_evidence",
    ]].reset_index(drop=True)
    values = dataset.values.loc[row_ids].reset_index(drop=True)
    return pd.concat([annotations, values], axis=1)
