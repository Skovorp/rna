#!/usr/bin/env python3
"""Crosswalk genes between two annotations of the same genome assembly.

The comparison is intentionally conservative: contigs are paired only when
their sequence length is unique in both references, and genes are paired only
when chromosome, start, end, and strand match exactly.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from dataclasses import dataclass
import gzip
from pathlib import Path
import re
from typing import TextIO


ATTRIBUTE_PATTERN = re.compile(r'(\S+)\s+"([^"]*)"')


@dataclass(frozen=True)
class Gene:
    contig: str
    start: int
    end: int
    strand: str
    gene_id: str
    gene_name: str


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def read_fai(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"Malformed FAI line in {path}: {line!r}")
            lengths[fields[0]] = int(fields[1])
    if not lengths:
        raise ValueError(f"No sequences found in {path}")
    return lengths


def read_genes(path: Path) -> list[Gene]:
    genes: list[Gene] = []
    seen_ids: set[str] = set()
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attributes = dict(ATTRIBUTE_PATTERN.findall(fields[8]))
            gene_id = attributes.get("gene_id", "").strip()
            if not gene_id:
                raise ValueError(f"Gene without gene_id at {path}:{line_number}")
            if gene_id in seen_ids:
                raise ValueError(f"Duplicate gene_id {gene_id!r} in {path}")
            seen_ids.add(gene_id)
            genes.append(
                Gene(
                    contig=fields[0],
                    start=int(fields[3]),
                    end=int(fields[4]),
                    strand=fields[6],
                    gene_id=gene_id,
                    gene_name=(
                        attributes.get("gene")
                        or attributes.get("gene_name")
                        or gene_id
                    ).strip(),
                )
            )
    if not genes:
        raise ValueError(f"No gene records found in {path}")
    return genes


def unique_length_contig_map(
    source_lengths: dict[str, int],
    target_lengths: dict[str, int],
) -> dict[str, str]:
    source_by_length: dict[int, list[str]] = defaultdict(list)
    target_by_length: dict[int, list[str]] = defaultdict(list)
    for contig, length in source_lengths.items():
        source_by_length[length].append(contig)
    for contig, length in target_lengths.items():
        target_by_length[length].append(contig)
    return {
        source_contigs[0]: target_by_length[length][0]
        for length, source_contigs in source_by_length.items()
        if len(source_contigs) == 1 and len(target_by_length.get(length, [])) == 1
    }


def write_crosswalk(
    output: Path,
    source_genes: list[Gene],
    target_genes: list[Gene],
    contig_map: dict[str, str],
) -> tuple[int, int]:
    target_index: dict[tuple[str, int, int, str], list[Gene]] = defaultdict(list)
    for gene in target_genes:
        target_index[(gene.contig, gene.start, gene.end, gene.strand)].append(gene)

    rows: list[dict[str, object]] = []
    ambiguous = 0
    for source in source_genes:
        target_contig = contig_map.get(source.contig)
        if target_contig is None:
            continue
        candidates = target_index.get(
            (target_contig, source.start, source.end, source.strand), []
        )
        if len(candidates) != 1:
            ambiguous += len(candidates) > 1
            continue
        target = candidates[0]
        rows.append(
            {
                "source_gene_id": source.gene_id,
                "source_gene_name": source.gene_name,
                "target_gene_id": target.gene_id,
                "target_gene_name": target.gene_name,
                "source_contig": source.contig,
                "target_contig": target.contig,
                "start": source.start,
                "end": source.end,
                "strand": source.strand,
                "match_method": "unique_contig_length_and_exact_gene_coordinates",
            }
        )

    target_counts = Counter(str(row["target_gene_id"]) for row in rows)
    duplicate_targets = {gene_id for gene_id, count in target_counts.items() if count > 1}
    if duplicate_targets:
        rows = [
            row for row in rows if str(row["target_gene_id"]) not in duplicate_targets
        ]

    output.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if output.suffix == ".gz" else open
    with opener(output, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows), ambiguous + len(duplicate_targets)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-gtf", type=Path, required=True)
    parser.add_argument("--source-fai", type=Path, required=True)
    parser.add_argument("--target-gtf", type=Path, required=True)
    parser.add_argument("--target-fai", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source_lengths = read_fai(args.source_fai)
    target_lengths = read_fai(args.target_fai)
    if Counter(source_lengths.values()) != Counter(target_lengths.values()):
        raise ValueError("Source and target references do not have matching sequence lengths")
    contig_map = unique_length_contig_map(source_lengths, target_lengths)
    source_genes = read_genes(args.source_gtf)
    target_genes = read_genes(args.target_gtf)
    matched, ambiguous = write_crosswalk(
        args.output, source_genes, target_genes, contig_map
    )
    print(
        f"Mapped {len(contig_map):,}/{len(source_lengths):,} contigs and "
        f"{matched:,}/{len(source_genes):,} source genes; "
        f"skipped {ambiguous:,} ambiguous exact-coordinate matches."
    )


if __name__ == "__main__":
    main()
