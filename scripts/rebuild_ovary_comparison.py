#!/usr/bin/env python3
"""Regenerate the ovary paper-vs-reprocessed comparison against current data.

`analysis/compare_elife_tpm.py` expects both matrices in the published
supplement's shape: `IDs`/`Symbols` columns followed by identically named
sample columns. The reprocessed STAR + Salmon matrix uses `gene_id`/`gene_name`
and nf-core sample names, so this adapts it before delegating to that script.

Run this whenever the reprocessed ovary matrix changes, or the comparison page
will keep describing a matrix the atlas no longer displays. Afterwards run
`scripts/theme_comparison_reports.py` and copy the full report into
`app/assets/ovary_comparison/`.

The app shows only `elife_ovary_tpm_full_report.html`. The generator also emits
a standalone zero-transition report, but that is a strict subset of the full
one (same figure, without the per-gene tables), so it is not bundled.
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PUBLISHED = ROOT / "expression" / "elife_80489_tpm.tsv.gz"
REPROCESSED = ROOT / "expression" / "ovary_star_salmon_gene_tpm.tsv.gz"
CROSSWALK = ROOT / "analysis" / "results" / "elife_tpm_comparison" / "annotation_crosswalk.tsv.gz"
OUTPUT_DIR = ROOT / "analysis" / "results" / "elife_tpm_comparison"
ADAPTED = ROOT / ".staging" / "ovary_reprocessed_as_published_shape.tsv.gz"


def published_sample_key(column: str) -> tuple[str, str]:
    """Fe_Ov_NBF_1_TPM -> ("NBF", "1")."""
    match = re.fullmatch(r"Fe_Ov_(.+)_([0-9]+)_TPM", column)
    if not match:
        raise ValueError(f"Unrecognized published sample column: {column}")
    return match.group(1), match.group(2)


def reprocessed_sample_key(column: str) -> tuple[str, str]:
    """Fe.Ov.12hBF.1_S10 -> ("12hBF", "1")."""
    match = re.fullmatch(r"Fe\.Ov\.(.+)\.([0-9]+)_S[0-9]+", column)
    if not match:
        raise ValueError(f"Unrecognized reprocessed sample column: {column}")
    condition, replicate = match.group(1), match.group(2)
    # The published supplement spells the two day-6 states without a dot.
    condition = condition.replace("6dBF.Retained", "6dBF_Retained")
    condition = condition.replace("6dBF.Laid", "6dBF_Laid")
    return condition, replicate


def main() -> None:
    published_columns = pd.read_csv(PUBLISHED, sep="\t", nrows=0).columns.tolist()[2:]
    reprocessed = pd.read_csv(REPROCESSED, sep="\t", low_memory=False)
    reprocessed_columns = reprocessed.columns.tolist()[2:]

    published_by_key = {published_sample_key(c): c for c in published_columns}
    reprocessed_by_key = {reprocessed_sample_key(c): c for c in reprocessed_columns}

    missing = sorted(set(published_by_key) - set(reprocessed_by_key))
    extra = sorted(set(reprocessed_by_key) - set(published_by_key))
    if missing or extra:
        raise SystemExit(
            "Sample sets do not correspond.\n"
            f"  only in published:   {missing}\n"
            f"  only in reprocessed: {extra}"
        )

    # Emit in the published column order so the comparison script's identity
    # check passes and samples line up positionally.
    adapted = pd.DataFrame(
        {
            "IDs": reprocessed["gene_id"],
            "Symbols": reprocessed["gene_name"].fillna(""),
        }
    )
    for key_column in published_columns:
        source = reprocessed_by_key[published_sample_key(key_column)]
        adapted[key_column] = reprocessed[source]

    ADAPTED.parent.mkdir(parents=True, exist_ok=True)
    adapted.to_csv(ADAPTED, sep="\t", index=False, compression="gzip")
    print(
        f"adapted reprocessed matrix: {adapted.shape[0]} genes x "
        f"{len(published_columns)} samples -> {ADAPTED}"
    )

    command = [
        sys.executable,
        str(ROOT / "analysis" / "compare_elife_tpm.py"),
        "--published", str(PUBLISHED),
        "--reanalysis", str(ADAPTED),
        "--crosswalk", str(CROSSWALK),
        "--output-dir", str(OUTPUT_DIR),
    ]
    print("running:", " ".join(command))
    subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
