#!/usr/bin/env python3
"""Build the app's ovary TPM matrix from per-sample Salmon gene quants."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
from pathlib import Path


REQUIRED_QUANT_COLUMNS = {"Name", "TPM"}
TPM_SUM_TOLERANCE = 0.5


def read_marker(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key] = value
    return values


def discover_quants(quant_root: Path) -> dict[str, Path]:
    quants: dict[str, Path] = {}
    for marker_path in sorted(quant_root.rglob("quant.complete.txt")):
        marker = read_marker(marker_path)
        sample_accession = marker.get("sample_accession", "")
        quant_path = marker_path.with_name("quant.genes.sf")
        if not sample_accession or not quant_path.is_file():
            continue
        if sample_accession in quants:
            raise ValueError(f"Duplicate completed quant for {sample_accession}")
        quants[sample_accession] = quant_path
    if not quants:
        raise ValueError(f"No completed gene quants found under {quant_root}")
    return quants


def read_samples(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"sample", "sample_accession"}
    if not rows or not required <= set(rows[0]):
        raise ValueError(f"{path} must contain {sorted(required)}")
    if len({row["sample"] for row in rows}) != len(rows):
        raise ValueError("Sample names must be unique")
    if len({row["sample_accession"] for row in rows}) != len(rows):
        raise ValueError("Sample accessions must be unique")
    return rows


def read_quant(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not REQUIRED_QUANT_COLUMNS <= set(reader.fieldnames or []):
            raise ValueError(f"{path} lacks Salmon Name/TPM columns")
        for row in reader:
            gene = row["Name"]
            if gene in values:
                raise ValueError(f"Duplicate gene {gene} in {path}")
            values[gene] = row["TPM"]
    tpm_sum = sum(float(value) for value in values.values())
    if abs(tpm_sum - 1_000_000.0) > TPM_SUM_TOLERANCE:
        raise ValueError(f"{path} TPM sum is {tpm_sum:.6f}, expected 1,000,000")
    return values


def write_matrix(
    output: Path,
    samples: list[dict[str, str]],
    values_by_sample: list[dict[str, str]],
) -> str:
    gene_sets = [set(values) for values in values_by_sample]
    if any(genes != gene_sets[0] for genes in gene_sets[1:]):
        raise ValueError("Salmon gene sets differ between biological samples")

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as raw_handle:
        with gzip.GzipFile(fileobj=raw_handle, mode="wb", filename="", mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as text_handle:
                writer = csv.writer(text_handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["IDs", "Symbols", *(row["sample"] for row in samples)])
                for gene in sorted(gene_sets[0]):
                    writer.writerow(
                        [gene, gene, *(values[gene] for values in values_by_sample)]
                    )
    return hashlib.sha256(output.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("quant_root", type=Path)
    parser.add_argument("sample_metadata", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    quants = discover_quants(args.quant_root)
    samples = read_samples(args.sample_metadata)
    expected_accessions = {row["sample_accession"] for row in samples}
    if set(quants) != expected_accessions:
        raise ValueError(
            "Quant/sample mismatch: "
            f"missing={sorted(expected_accessions - set(quants))}, "
            f"unexpected={sorted(set(quants) - expected_accessions)}"
        )

    for row in samples:
        marker = read_marker(
            quants[row["sample_accession"]].with_name("quant.complete.txt")
        )
        expected_alias = row.get("sample_alias", "")
        if expected_alias and marker.get("sample_alias") != expected_alias:
            raise ValueError(
                f"Alias mismatch for {row['sample_accession']}: "
                f"metadata={expected_alias}, marker={marker.get('sample_alias', '')}"
            )

    values_by_sample = [
        read_quant(quants[row["sample_accession"]]) for row in samples
    ]
    sha256 = write_matrix(args.output, samples, values_by_sample)
    print(
        f"Wrote {len(values_by_sample[0]):,} genes x {len(samples)} samples to "
        f"{args.output} (sha256 {sha256})"
    )


if __name__ == "__main__":
    main()
