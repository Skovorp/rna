#!/usr/bin/env python3
"""Build the Matthews 2016 paper-vs-reprocessed tissue-atlas comparison."""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PUBLISHED = ROOT / "expression" / "neurotranscriptome_2016_aaegl_ru_tpm.tsv.gz"
REPROCESSED = ROOT / "expression" / "atlas_star_salmon_gene_tpm.tsv.gz"
SAMPLE_METADATA = ROOT / "expression" / "neurotranscriptome_2016_samples.tsv"
OUTPUT_DIR = ROOT / "analysis" / "results" / "atlas_tpm_comparison"
STAGING_DIR = ROOT / ".staging" / "atlas_comparison"
ASSET_DIR = ROOT / "app" / "assets" / "atlas_comparison"
REPORT_FILENAME = "matthews_2016_atlas_tpm_full_report.html"


def normalize_sample(value: str) -> str:
    """Match case-only differences and nf-core's duplicate-name suffixes."""
    return re.sub(r"\.[0-9]+$", "", value).casefold()


def joined_symbols(values: pd.Series) -> str:
    symbols = [
        str(value).strip()
        for value in values
        if pd.notna(value) and str(value).strip()
    ]
    return " / ".join(dict.fromkeys(symbols))


def adapt_published() -> tuple[pd.DataFrame, list[str]]:
    paper = pd.read_csv(PUBLISHED, sep="\t", low_memory=False)
    identifier = "Vectorbase Identifier"
    sample_columns = paper.columns.tolist()[3:]
    if not sample_columns:
        raise ValueError("Published tissue-atlas matrix has no sample columns")
    if (
        paper[identifier].isna().any()
        or paper[identifier].astype(str).str.strip().eq("").any()
    ):
        raise ValueError("Published tissue-atlas matrix has missing identifiers")

    # The paper matrix contains 22 historical identifiers twice. Its rows are
    # transcript-model features while our matrix is gene-level, so sum TPM for
    # the repeated stable identifier before making a one-to-one comparison.
    numeric = paper[sample_columns].apply(pd.to_numeric, errors="raise")
    numeric.insert(0, identifier, paper[identifier].astype(str))
    collapsed_values = numeric.groupby(identifier, sort=False, as_index=False).sum()
    collapsed_symbols = (
        paper.groupby(identifier, sort=False)["Display name"]
        .agg(joined_symbols)
        .rename("Symbols")
        .reset_index()
    )
    adapted = collapsed_symbols.merge(
        collapsed_values, on=identifier, how="inner", validate="one_to_one"
    ).rename(columns={identifier: "IDs"})
    adapted = adapted[["IDs", "Symbols", *sample_columns]]
    return adapted, sample_columns


def adapt_reprocessed(published_samples: list[str]) -> pd.DataFrame:
    reprocessed = pd.read_csv(REPROCESSED, sep="\t", low_memory=False)
    if reprocessed["gene_id"].duplicated().any():
        raise ValueError("Reprocessed tissue-atlas matrix has duplicate gene IDs")

    reprocessed_samples = reprocessed.columns.tolist()[2:]
    by_normalized = {normalize_sample(sample): sample for sample in reprocessed_samples}
    if len(by_normalized) != len(reprocessed_samples):
        raise ValueError("Reprocessed tissue-atlas sample names are ambiguous")
    missing = [
        sample
        for sample in published_samples
        if normalize_sample(sample) not in by_normalized
    ]
    if missing:
        raise ValueError(f"Published samples missing from reprocessing: {missing}")

    return pd.DataFrame(
        {
            "IDs": reprocessed["gene_id"].astype(str),
            "Symbols": reprocessed["gene_name"].fillna("").astype(str),
            **{
                sample: reprocessed[by_normalized[normalize_sample(sample)]].to_numpy()
                for sample in published_samples
            },
        }
    )


def adapt_metadata(published_samples: list[str]) -> pd.DataFrame:
    metadata = pd.read_csv(SAMPLE_METADATA, sep="\t", dtype=str).fillna("")
    metadata = metadata.set_index("library_id", drop=False)
    missing = sorted(set(published_samples) - set(metadata.index))
    if missing:
        raise ValueError(f"Published sample metadata is missing: {missing}")
    metadata = metadata.loc[published_samples].copy()
    metadata["sample"] = metadata["library_id"]
    metadata["tissue"] = metadata["tissue"].str.replace(
        "abdominaltip", "abdominal tip", regex=False
    )
    metadata["reproductive_state"] = metadata["tissue"].str.title()
    metadata["replicate"] = metadata["sample"].str.extract(r"_([0-9]+)$")[0]
    return metadata[["sample", "reproductive_state", "replicate"]]


def main() -> None:
    published, published_samples = adapt_published()
    reprocessed = adapt_reprocessed(published_samples)
    metadata = adapt_metadata(published_samples)

    STAGING_DIR.mkdir(parents=True, exist_ok=True)
    published_path = STAGING_DIR / "published_gene_tpm.tsv.gz"
    reprocessed_path = STAGING_DIR / "reprocessed_gene_tpm.tsv.gz"
    metadata_path = STAGING_DIR / "samples.tsv"
    crosswalk_path = STAGING_DIR / "empty_crosswalk.tsv"
    published.to_csv(published_path, sep="\t", index=False, compression="gzip")
    reprocessed.to_csv(
        reprocessed_path, sep="\t", index=False, compression="gzip"
    )
    metadata.to_csv(metadata_path, sep="\t", index=False)
    pd.DataFrame(
        columns=["source_gene_id", "source_gene_name", "target_gene_id"]
    ).to_csv(crosswalk_path, sep="\t", index=False)

    print(
        f"prepared {len(published):,} unique published genes and "
        f"{len(published_samples)} matched samples"
    )
    command = [
        sys.executable,
        str(ROOT / "analysis" / "compare_elife_tpm.py"),
        "--published",
        str(published_path),
        "--reanalysis",
        str(reprocessed_path),
        "--metadata",
        str(metadata_path),
        "--crosswalk",
        str(crosswalk_path),
        "--output-dir",
        str(OUTPUT_DIR),
        "--group-label",
        "Tissue",
        "--report-title",
        "Published versus reanalysed tissue-atlas TPM",
        "--report-subtitle",
        (
            "Matthews et al. BMC Genomics 2016, 122 matched biological samples, "
            "log expression is <code>log2(TPM + 1)</code>."
        ),
        "--report-filename",
        REPORT_FILENAME,
    ]
    subprocess.run(command, check=True)

    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    shutil.copy2(OUTPUT_DIR / "figures.json", ASSET_DIR / "figures.json")
    shutil.copy2(OUTPUT_DIR / REPORT_FILENAME, ASSET_DIR / REPORT_FILENAME)
    subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "theme_comparison_reports.py")],
        check=True,
    )
    print(f"prepared live comparison assets in {ASSET_DIR}")


if __name__ == "__main__":
    main()
