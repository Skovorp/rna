"""Load, harmonize, search, and summarize the bundled TPM matrices."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path
import pickle
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


DATASET_ORDER = (
    "ovary_paper",
    "elife",
    "neuro_ru",
    "midgut",
    "crop",
    "neuro_legacy",
)

STAR_SALMON_CONDITION_LABELS = {
    "NBF": "Non-blood-fed",
    "3hBF": "3 hours post-blood-meal",
    "6hBF": "6 hours post-blood-meal",
    "12hBF": "12 hours post-blood-meal",
    "24hBF": "24 hours post-blood-meal",
    "48hBF": "48 hours post-blood-meal",
    "72hBF": "72 hours post-blood-meal",
    "96hBF": "96 hours post-blood-meal",
    "6dBF.Retained": "6 days post-blood-meal (eggs retained)",
    "6dBF.Laid": "6 days post-blood-meal (eggs laid)",
    "13dBF": "13 days post-blood-meal",
    "Ma.Mg": "Male, non-blood-fed",
}

STAR_SALMON_CONDITION_SEQUENCE = tuple(STAR_SALMON_CONDITION_LABELS)

# Elapsed time since the blood meal, independent of sex. Non-blood-fed samples
# and the male midgut have no post-blood-meal time.
STAR_SALMON_TIMEPOINT_LABELS = {
    "NBF": "not_applicable",
    "3hBF": "3h",
    "6hBF": "6h",
    "12hBF": "12h",
    "24hBF": "24h",
    "48hBF": "48h",
    "72hBF": "72h",
    "96hBF": "96h",
    "6dBF.Retained": "6d",
    "6dBF.Laid": "6d",
    "13dBF": "13d",
    "Ma.Mg": "not_applicable",
}

CONDITION_LABELS = {
    "BF": "Blood-fed",
    "O": "Gravid / oviposition-stage",
    "SF": "Non-blood-fed / sugar-fed",
    "male": "Male",
}

ELIFE_CONDITION_IDS = {
    "Non-blood-fed": "female_NBF",
    "3 hours post-blood-meal": "female_3hBF",
    "6 hours post-blood-meal": "female_6hBF",
    "12 hours post-blood-meal": "female_12hBF",
    "24 hours (1 day) post-blood-meal": "female_24hBF",
    "48 hours (2 days) post-blood-meal": "female_48hBF",
    "72 hours (3 days) post-blood-meal": "female_72hBF",
    "96 hours (4 days) post-blood-meal": "female_96hBF",
    "6 days post-blood-meal (eggs retained)": "female_6dBF_retained",
    "6 days post-blood-meal (eggs laid <5 hours prior)": "female_6dBF_laid",
    "13 days post-blood-meal (eggs laid >1 week prior)": "female_13dBF",
}

PAPER_FAMILY_LABELS = {
    "IR": "Ionotropic receptors (IR)",
    "OR": "Odorant receptors (OR)",
    "GR": "Gustatory receptors (GR)",
    "OBP": "Odorant-binding proteins (OBP)",
}

ORCO_ALIASES = ("Orco", "AaegOr7", "Or7", "AAEL005776")

# Bump when a change alters derived dataset content without touching this file.
_CACHE_VERSION = 1
_CACHE_DIR = ".cache"

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


_NON_ALPHANUMERIC = re.compile(r"[^a-z0-9]")


def _clean_series(values: pd.Series) -> np.ndarray:
    """Vectorized _clean over a whole column."""
    return values.fillna("").astype(str).str.strip().to_numpy()


def normalize_alias(value: object) -> str:
    """Normalize an identifier for forgiving exact matching."""
    text = _clean(value).casefold()
    return _NON_ALPHANUMERIC.sub("", text)


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


# Compiled once: these were rebuilt from an f-string on every call, which cost
# ~440k regex-cache lookups per startup.
_NUMBERED_SYMBOL = r"\d+(?:[A-Z]+\d*)?"
_FAMILY_PATTERNS = (
    (re.compile(rf"OR{_NUMBERED_SYMBOL}"), "Odorant receptors (OR)"),
    (re.compile(rf"IR{_NUMBERED_SYMBOL}"), "Ionotropic receptors (IR)"),
    (re.compile(rf"GR{_NUMBERED_SYMBOL}"), "Gustatory receptors (GR)"),
    (re.compile(rf"OBP{_NUMBERED_SYMBOL}"), "Odorant-binding proteins (OBP)"),
)


def classify_family(symbol: object) -> str:
    normalized = normalize_alias(symbol).upper()
    if normalized == "ORCO":
        return "Odorant receptors (OR)"
    for pattern, family in _FAMILY_PATTERNS:
        if pattern.fullmatch(normalized):
            return family
    return "Other"


def _finalize_genes(
    annotations: pd.DataFrame,
    stable_id: pd.Series,
    internal_id: pd.Series,
    raw_symbol: pd.Series,
) -> pd.DataFrame:
    genes = annotations.copy().reset_index(drop=True)
    genes.insert(0, "row_id", np.arange(len(genes), dtype=int))
    genes["stable_id"] = _clean_series(stable_id)
    genes["internal_id"] = _clean_series(internal_id)
    genes["raw_symbol"] = _clean_series(raw_symbol)
    genes["canonical_symbol"] = [
        canonical_symbol(symbol, stable)
        for symbol, stable in zip(genes["raw_symbol"], genes["stable_id"])
    ]
    genes["display_name"] = genes["canonical_symbol"].where(
        genes["canonical_symbol"].ne(""), genes["stable_id"]
    )
    # Symbols repeat heavily across a matrix, so classify each distinct one once.
    family_by_symbol = {
        symbol: classify_family(symbol)
        for symbol in genes["canonical_symbol"].unique()
    }
    genes["family"] = genes["canonical_symbol"].map(family_by_symbol)
    for column in PAPER_ANNOTATION_COLUMNS:
        genes[column] = ""

    # Built by zipping raw columns rather than genes.apply(axis=1): the latter
    # materializes a Series per gene and dominated startup (~560k Series
    # lookups across the bundled matrices).
    alias_rows: list[tuple[str, ...]] = []
    search_text: list[str] = []
    search_normalized: list[str] = []
    normalized_cache: dict[str, str] = {}
    for stable, internal, raw, canonical in zip(
        genes["stable_id"],
        genes["internal_id"],
        genes["raw_symbol"],
        genes["canonical_symbol"],
    ):
        values = {stable, internal, raw, canonical}
        if canonical.casefold() == "orco":
            values.update(ORCO_ALIASES)
        row = tuple(sorted(value for value in values if value))
        alias_rows.append(row)
        search_text.append(" | ".join(row))
        normalized = []
        for value in row:
            cached = normalized_cache.get(value)
            if cached is None:
                cached = normalize_alias(value)
                normalized_cache[value] = cached
            normalized.append(cached)
        search_normalized.append(" | ".join(normalized))

    genes["aliases"] = alias_rows
    genes["search_text"] = search_text
    genes["search_normalized"] = search_normalized
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
        samples["tissue"].str.title() + ", " + samples["condition_label"]
    )
    samples["reproductive_state"] = samples["condition_label"]
    samples = samples.set_index("sample", drop=False).reindex(list(sample_columns))
    return samples


def _load_elife_samples(path: Path, sample_columns: Iterable[str]) -> pd.DataFrame:
    samples = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    samples["sample"] = samples["sample"].astype(str)
    samples["condition"] = samples["reproductive_state"].map(ELIFE_CONDITION_IDS)
    if samples["condition"].isna().any():
        unknown = sorted(
            samples.loc[
                samples["condition"].isna(), "reproductive_state"
            ].unique()
        )
        raise ValueError(f"Unknown eLife reproductive states: {unknown}")
    samples["condition_label"] = samples["reproductive_state"]
    samples["tissue"] = "ovary"
    samples["sex"] = "female"
    samples["tissue_condition"] = "Ovary, " + samples["reproductive_state"]
    samples = samples.set_index("sample", drop=False).reindex(list(sample_columns))
    return samples


def _load_midgut_samples(path: Path, sample_columns: Iterable[str]) -> pd.DataFrame:
    samples = pd.read_csv(path, dtype=str).fillna("")
    condition_labels = {
        "male_midgut": "Male, non-blood-fed",
        "female_NBF": "Female, non-blood-fed",
        "female_3hBF": "Female, 3 h post-blood-meal",
        "female_6hBF": "Female, 6 h post-blood-meal",
        "female_12hBF": "Female, 12 h post-blood-meal",
        "female_24hBF": "Female, 24 h post-blood-meal",
        "female_48hBF": "Female, 48 h post-blood-meal",
        "female_72hBF": "Female, 72 h post-blood-meal",
    }
    samples["condition_label"] = samples["condition"].map(condition_labels).fillna(
        samples["condition"]
    )
    samples["reproductive_state"] = samples["condition_label"]
    samples["tissue"] = "midgut"
    samples["tissue_condition"] = "Midgut, " + samples["condition_label"]
    condition_order = {
        "female_NBF": 0,
        "female_3hBF": 1,
        "female_6hBF": 2,
        "female_12hBF": 3,
        "female_24hBF": 4,
        "female_48hBF": 5,
        "female_72hBF": 6,
        "male_midgut": 7,
    }
    samples = samples[samples["sample"].isin(sample_columns)].copy()
    samples["_condition_order"] = samples["condition"].map(condition_order).fillna(99)
    samples["_replicate_order"] = pd.to_numeric(samples["replicate"], errors="coerce")
    samples = samples.sort_values(
        ["_condition_order", "_replicate_order", "sample"], kind="stable"
    ).drop(columns=["_condition_order", "_replicate_order"])
    samples = samples.set_index("sample", drop=False)
    return samples


def star_salmon_condition(sample: str) -> str:
    """Fe.Mg.12hBF.1_S13 / Fe_Crop_NBF_1_S1 -> 12hBF / NBF; Ma.Mg.1_S1 -> Ma.Mg."""
    condition = re.sub(r"^Fe[._][A-Za-z]+[._]", "", sample)
    return re.sub(r"[._][0-9]+_S[0-9]+$", "", condition)


def _load_star_salmon_samples(
    sample_columns: Iterable[str], tissue: str
) -> pd.DataFrame:
    samples = pd.DataFrame({"sample": list(sample_columns)})
    samples["condition"] = samples["sample"].map(star_salmon_condition)
    unknown = sorted(
        set(samples["condition"]) - set(STAR_SALMON_CONDITION_LABELS)
    )
    if unknown:
        raise ValueError(f"Unknown star_salmon conditions: {unknown}")
    samples["condition_label"] = samples["condition"].map(
        STAR_SALMON_CONDITION_LABELS
    )
    samples["reproductive_state"] = samples["condition_label"]
    samples["sex"] = np.where(
        samples["sample"].str.startswith("Ma"), "male", "female"
    )
    samples["tissue"] = tissue
    samples["tissue_condition"] = (
        tissue.title() + ", " + samples["condition_label"]
    )
    samples["replicate"] = samples["sample"].str.extract(
        r"[._]([0-9]+)_S[0-9]+$"
    )[0]
    samples["timepoint"] = samples["condition"].map(
        STAR_SALMON_TIMEPOINT_LABELS
    ).fillna("not_applicable")
    samples["_condition_order"] = samples["condition"].map(
        {value: index for index, value in enumerate(STAR_SALMON_CONDITION_SEQUENCE)}
    )
    samples["_replicate_order"] = pd.to_numeric(
        samples["replicate"], errors="coerce"
    )
    samples = samples.sort_values(
        ["_condition_order", "_replicate_order", "sample"], kind="stable"
    ).drop(columns=["_condition_order", "_replicate_order"])
    return samples.set_index("sample", drop=False)


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
    *,
    label: str | None = None,
    paper: str = "Local nf-core/rnaseq import",
    annotation_version: str = "Identifiers from imported matrix",
    samples: pd.DataFrame | None = None,
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
    matrix_symbols = (
        frame["gene_name"].map(_clean)
        if "gene_name" in frame.columns
        else pd.Series([""] * len(frame))
    )
    symbols = stable_ids.map(crosswalk).fillna("")
    symbols = symbols.where(symbols.ne(""), matrix_symbols)
    annotation_columns = [gene_column]
    if "gene_name" in frame.columns and "gene_name" != gene_column:
        annotation_columns.append("gene_name")
    annotations = frame[annotation_columns].fillna("")
    genes = _finalize_genes(
        annotations,
        stable_ids,
        pd.Series([""] * len(frame)),
        symbols,
    )
    values = frame.drop(columns=annotation_columns).apply(pd.to_numeric, errors="coerce")
    if values.isna().all(axis=None):
        raise ValueError(f"No numeric TPM values found in {path}")
    values = values.fillna(0.0)
    values.index = np.arange(len(values), dtype=int)
    if label is None:
        label = path.name
        for suffix in (".gz", ".tsv"):
            if label.endswith(suffix):
                label = label[: -len(suffix)]
        label = f"nf-core {label}"
    if samples is None:
        samples = _generic_samples(values.columns)
    else:
        matrix_samples = set(values.columns)
        metadata_samples = set(samples.index)
        if matrix_samples != metadata_samples:
            missing_metadata = sorted(matrix_samples - metadata_samples)
            missing_values = sorted(metadata_samples - matrix_samples)
            raise ValueError(
                "TPM/sample metadata mismatch: "
                f"missing metadata={missing_metadata}, missing TPM={missing_values}"
            )
        values = values.loc[:, samples.index]
    return ExpressionDataset(
        key=key,
        label=label,
        paper=paper,
        annotation_version=annotation_version,
        genes=genes,
        values=values,
        samples=samples,
    )


def _cache_signature(expression_dir: Path) -> str:
    """Identify the inputs and the code that derive the dataset objects."""
    parts = [f"v{_CACHE_VERSION}", str(Path(__file__).stat().st_mtime_ns)]
    for path in sorted(expression_dir.rglob("*")):
        if path.is_file() and not path.is_relative_to(expression_dir / _CACHE_DIR):
            stat = path.stat()
            parts.append(f"{path.name}:{stat.st_size}:{stat.st_mtime_ns}")
    return hashlib.sha256("|".join(parts).encode()).hexdigest()[:16]


def load_datasets(
    expression_dir: Path | str, use_cache: bool = True
) -> dict[str, ExpressionDataset]:
    """Build the dataset objects, reusing a derived cache when inputs are unchanged.

    Parsing and harmonizing the bundled matrices takes a couple of seconds, which
    is paid on every server start and by every test. The cache is keyed by the
    size and mtime of every source file plus this module's own mtime, so editing
    either invalidates it automatically.
    """
    expression_dir = Path(expression_dir)
    cache_path = None
    if use_cache:
        cache_path = (
            expression_dir / _CACHE_DIR / f"datasets-{_cache_signature(expression_dir)}.pkl"
        )
        try:
            with cache_path.open("rb") as handle:
                return pickle.load(handle)
        except (OSError, EOFError, pickle.UnpicklingError, AttributeError):
            pass

    datasets = _build_datasets(expression_dir)

    if cache_path is not None:
        try:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            # Write-then-rename so a concurrent reader never sees a partial file.
            temporary = cache_path.with_suffix(".tmp")
            with temporary.open("wb") as handle:
                pickle.dump(datasets, handle, protocol=pickle.HIGHEST_PROTOCOL)
            temporary.replace(cache_path)
            for stale in cache_path.parent.glob("datasets-*.pkl"):
                if stale != cache_path:
                    stale.unlink(missing_ok=True)
        except OSError:
            pass  # A read-only bundle just pays the rebuild cost.
    return datasets


def _build_datasets(expression_dir: Path) -> dict[str, ExpressionDataset]:
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

    paper_annotations, paper_values = _read_matrix(
        expression_dir / "elife_80489_tpm.tsv.gz",
        ["IDs", "Symbols"],
    )
    paper_genes = _finalize_genes(
        paper_annotations,
        paper_annotations["IDs"],
        pd.Series([""] * len(paper_annotations)),
        paper_annotations["Symbols"],
    )
    paper_samples = _load_elife_samples(
        expression_dir / "elife_80489_samples.tsv", paper_values.columns
    )

    def star_salmon_dataset(
        file_name: str, key: str, tissue: str, label: str, paper: str
    ) -> ExpressionDataset:
        path = expression_dir / file_name
        columns = pd.read_csv(path, sep="\t", compression="gzip", nrows=0).columns
        sample_columns = [
            column for column in columns if column not in {"gene_id", "gene_name"}
        ]
        samples = _load_star_salmon_samples(sample_columns, tissue)
        return load_nfcore_dataset(
            path,
            key,
            crosswalk,
            label=label,
            paper=paper,
            annotation_version="AaegL5, VectorBase 58 + Jové et al. 2019",
            samples=samples,
        )

    ovary = star_salmon_dataset(
        "ovary_star_salmon_gene_tpm.tsv.gz",
        "elife",
        "ovary",
        "Ovary (reprocessed), blood-meal time course",
        "Venkataraman et al. raw reads, our nf-core reprocessing",
    )
    midgut = star_salmon_dataset(
        "midgut_star_salmon_gene_tpm.tsv.gz",
        "midgut",
        "midgut",
        "Midgut (reprocessed), blood-meal time course",
        "Nadav Shai, Vosshall lab midgut RNA-seq",
    )
    crop = star_salmon_dataset(
        "crop_star_salmon_gene_tpm.tsv.gz",
        "crop",
        "crop",
        "Crop (reprocessed), non-blood-fed",
        "Vosshall lab crop RNA-seq",
    )

    datasets = {
        "ovary_paper": ExpressionDataset(
            key="ovary_paper",
            label="Ovary (paper), published TPM",
            paper="Venkataraman et al., eLife 2023",
            annotation_version="Published TPM supplement",
            genes=paper_genes,
            values=paper_values,
            samples=paper_samples,
        ),
        "elife": ovary,
        "neuro_ru": ExpressionDataset(
            key="neuro_ru",
            label="Atlas (paper), neurotranscriptome AaegL.RU",
            paper="Matthews et al., BMC Genomics 2016",
            annotation_version="AaegL.RU (recommended)",
            genes=ru_genes,
            values=ru_values,
            samples=bmc_samples,
        ),
        "neuro_legacy": ExpressionDataset(
            key="neuro_legacy",
            label="Atlas (paper), neurotranscriptome legacy AaegL3.3",
            paper="Matthews et al., BMC Genomics 2016",
            annotation_version="AaegL3.3 (compatibility)",
            genes=legacy_genes,
            values=legacy_values,
            samples=legacy_samples,
        ),
        "midgut": midgut,
        "crop": crop,
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
