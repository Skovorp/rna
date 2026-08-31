"""Resolve local gene aliases to embeddable UCSC Cell Browser views."""

from __future__ import annotations

from dataclasses import dataclass
import gzip
import json
import os
from pathlib import Path
import re
from urllib.parse import urlencode


UCSC_CELL_BROWSER_DEFAULT = "https://cells.ucsc.edu/"


def _cell_browser_base() -> str:
    """Use the public UCSC host unless deployment config selects our proxy."""
    configured = os.environ.get(
        "UCSC_CELL_BROWSER_BASE", UCSC_CELL_BROWSER_DEFAULT
    ).strip()
    return f"{configured.rstrip('/')}/"


def _normalize(value: object) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).strip().casefold())


@dataclass(frozen=True)
class AtlasDataset:
    name: str
    label: str
    parent_label: str
    sample_count: int
    aliases: dict[str, tuple[str, ...]]
    quick_genes: tuple[tuple[str, str], ...]
    categorical_fields: tuple[tuple[str, str], ...]
    default_metadata_field: str

    @property
    def display_label(self) -> str:
        label = self.label
        if self.parent_label:
            label = f"{self.parent_label}: {label}"
        return f"{label}, {self.sample_count:,} nuclei"


@dataclass(frozen=True)
class AtlasGeneMatch:
    dataset: AtlasDataset
    gene_query: str


def load_manifest(path: Path | str) -> list[AtlasDataset]:
    """Load the checked-in UCSC dataset/gene index."""
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        payload = json.load(handle)

    datasets: list[AtlasDataset] = []
    for item in payload["datasets"]:
        alias_queries: dict[str, set[str]] = {}
        for indexed_name in item["genes"]:
            parts = indexed_name.split("|")
            gene_query = parts[0]
            for alias in (indexed_name, *parts):
                normalized = _normalize(alias)
                if normalized:
                    alias_queries.setdefault(normalized, set()).add(gene_query)
        aliases = {
            alias: tuple(sorted(queries))
            for alias, queries in alias_queries.items()
        }
        datasets.append(
            AtlasDataset(
                name=item["name"],
                label=item["label"],
                parent_label=item.get("parent_label", ""),
                sample_count=int(item["sample_count"]),
                aliases=aliases,
                quick_genes=tuple(
                    (str(quick_gene["gene"]), str(quick_gene.get("label", "")))
                    for quick_gene in item.get("quick_genes", [])
                ),
                categorical_fields=tuple(
                    (str(field["name"]), str(field.get("label", field["name"])))
                    for field in item.get("categorical_fields", [])
                ),
                default_metadata_field=str(item.get("default_metadata_field", "")),
            )
        )
    return datasets


def find_gene_matches(
    datasets: list[AtlasDataset], aliases: list[str] | tuple[str, ...]
) -> list[AtlasGeneMatch]:
    """Return atlas views accepting any exact local alias for a gene."""
    normalized_aliases = list(
        dict.fromkeys(_normalize(alias) for alias in aliases if _normalize(alias))
    )
    matches: list[AtlasGeneMatch] = []
    for dataset in datasets:
        for alias in normalized_aliases:
            gene_queries = dataset.aliases.get(alias, ())
            if len(gene_queries) == 1:
                matches.append(AtlasGeneMatch(dataset, gene_queries[0]))
                break
    return matches


def cell_browser_url(dataset_name: str, gene_query: str) -> str:
    """Build the shareable UCSC URL used by both links and iframes."""
    query = urlencode(
        {
            "ds": dataset_name.replace("/", " "),
            "gene": gene_query,
        }
    )
    return f"{_cell_browser_base()}?{query}"


def cell_browser_metadata_url(dataset_name: str, metadata_field: str) -> str:
    """Build a UCSC UMAP URL colored by one metadata field."""
    query = urlencode(
        {
            "ds": dataset_name.replace("/", " "),
            "meta": metadata_field,
        }
    )
    return f"{_cell_browser_base()}?{query}"


def cell_browser_expression_url(
    dataset_name: str,
    gene_queries: list[str] | tuple[str, ...],
    metadata_field: str,
    context_gene: str | None = None,
) -> str:
    """Build a UCSC multi-gene dot-plot URL for one atlas view."""
    parameters = {
        "ds": dataset_name.replace("/", " "),
    }
    if context_gene:
        parameters["gene"] = context_gene
    parameters.update(
        {
            "exprGene": " ".join(gene_queries),
            "exprMeta": metadata_field,
        }
    )
    return f"{_cell_browser_base()}?{urlencode(parameters)}"
