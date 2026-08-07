#!/usr/bin/env python3
"""Snapshot UCSC Mosquito Cell Atlas datasets and accepted gene names."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import io
import json
from pathlib import Path
import ssl
from urllib.request import urlopen


BASE_URL = "https://cells.ucsc.edu/"
ROOT_DATASET = "mosquito"
DEFAULT_OUTPUT = (
    Path(__file__).resolve().parents[1]
    / "expression"
    / "ucsc_mosquito_cell_atlas_genes.json.gz"
)


def ssl_context() -> ssl.SSLContext:
    context = ssl.create_default_context()
    try:
        import certifi
    except ImportError:
        return context
    context.load_verify_locations(certifi.where())
    return context


def fetch_json(relative_url: str, context: ssl.SSLContext) -> object:
    with urlopen(BASE_URL + relative_url, timeout=120, context=context) as response:
        return json.load(response)


def leaf_datasets(
    dataset_name: str,
    context: ssl.SSLContext,
    parent_label: str = "",
) -> list[dict[str, object]]:
    payload = fetch_json(f"{dataset_name}/dataset.json", context)
    children = payload.get("datasets", [])
    if not children:
        return [
            {
                "name": payload["name"],
                "label": payload.get("shortLabel", payload["name"]),
                "parent_label": parent_label,
                "sample_count": payload["sampleCount"],
            }
        ]

    leaves: list[dict[str, object]] = []
    for child in children:
        if child.get("isCollection"):
            leaves.extend(
                leaf_datasets(
                    child["name"],
                    context,
                    parent_label=child["shortLabel"],
                )
            )
        else:
            leaves.append(
                {
                    "name": child["name"],
                    "label": child["shortLabel"],
                    "parent_label": parent_label,
                    "sample_count": child["sampleCount"],
                }
            )
    return leaves


def build_manifest() -> dict[str, object]:
    context = ssl_context()
    datasets = leaf_datasets(ROOT_DATASET, context)

    def add_genes(dataset: dict[str, object]) -> dict[str, object]:
        genes = fetch_json(f"{dataset['name']}/exprMatrix.json", context)
        return {**dataset, "genes": sorted(genes)}

    with ThreadPoolExecutor(max_workers=8) as pool:
        indexed = list(pool.map(add_genes, datasets))

    indexed.sort(key=lambda item: (item["name"] != "mosquito/all", item["name"]))
    return {
        "source": f"{BASE_URL}{ROOT_DATASET}/dataset.json",
        "datasets": indexed,
    }


def write_manifest(payload: dict[str, object], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    encoded = json.dumps(payload, separators=(",", ":"), ensure_ascii=False).encode()
    buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode="wb", mtime=0) as archive:
        archive.write(encoded)
    output.write_bytes(buffer.getvalue())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = build_manifest()
    write_manifest(payload, args.output)
    gene_counts = [len(dataset["genes"]) for dataset in payload["datasets"]]
    print(
        f"wrote {args.output}: {len(gene_counts)} datasets, "
        f"{min(gene_counts):,}–{max(gene_counts):,} genes each"
    )


if __name__ == "__main__":
    main()
