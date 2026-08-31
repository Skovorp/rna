#!/usr/bin/env python3
"""Extract a unique transcript-to-gene mapping from an NCBI GTF file."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


ATTRIBUTE = re.compile(r'([A-Za-z0-9_]+) "([^"]*)"')


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("gtf", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    mappings: dict[str, str] = {}
    with args.gtf.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed GTF line {line_number}: expected 9 fields")
            attributes = dict(ATTRIBUTE.findall(fields[8]))
            transcript_id = attributes.get("transcript_id")
            gene_id = attributes.get("gene_id")
            if not transcript_id or not gene_id:
                continue
            previous = mappings.setdefault(transcript_id, gene_id)
            if previous != gene_id:
                raise ValueError(
                    f"Transcript {transcript_id} maps to both {previous} and {gene_id}"
                )

    if not mappings:
        raise ValueError("No transcript_id/gene_id pairs found")

    with args.output.open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["transcript_id", "gene_id"])
        writer.writerows(sorted(mappings.items()))

    print(f"Wrote {len(mappings):,} transcript-to-gene mappings to {args.output}")


if __name__ == "__main__":
    main()
