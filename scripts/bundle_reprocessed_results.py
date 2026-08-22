#!/usr/bin/env python3
"""Bundle the identically-processed nf-core star_salmon reruns into expression/.

Input: .staging/<dataset>/results/{star_salmon,de_pairwise} pulled from the Pi
rsync landing dir. Output: gzipped TPM matrices and DESeq2 contrast bundles in
the repository layout the app loads.
"""
from __future__ import annotations

import gzip
import re
import shutil
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
STAGING = ROOT / ".staging"
EXPRESSION = ROOT / "expression"

DATASETS = ("ovary", "midgut", "crop", "yedlin")


def condition_from_sample(name: str) -> str:
    """Fe.Mg.12hBF.1_S13 / Fe_Crop_NBF_1_S1 -> 12hBF / NBF; Ma.Mg.1_S1 -> Ma.Mg."""
    cond = re.sub(r"^Fe[._][A-Za-z]+[._]", "", name)
    cond = re.sub(r"[._][0-9]+_S[0-9]+$", "", cond)
    return cond


def bundle_tpm(dataset: str) -> list[str]:
    src = STAGING / dataset / "results" / "star_salmon" / "salmon.merged.gene_tpm.tsv"
    out = EXPRESSION / f"{dataset}_star_salmon_gene_tpm.tsv.gz"
    frame = pd.read_csv(src, sep="\t")
    assert list(frame.columns[:2]) == ["gene_id", "gene_name"], frame.columns[:5]
    with gzip.open(out, "wt") as handle:
        frame.to_csv(handle, sep="\t", index=False)
    samples = [c for c in frame.columns if c not in ("gene_id", "gene_name")]
    print(f"{dataset}: {len(frame)} genes x {len(samples)} samples -> {out.name}")
    return samples


def bundle_de(dataset: str) -> None:
    src_dir = STAGING / dataset / "results" / "de_pairwise"
    if not src_dir.is_dir():
        print(f"{dataset}: no de_pairwise (single condition), skipping DE")
        return
    out_dir = EXPRESSION / f"{dataset}_deseq2"
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)
    rows = []
    for path in sorted(src_dir.glob("DE_*.tsv")):
        match = re.fullmatch(r"DE_(.+)_vs_(.+)\.tsv", path.name)
        target, reference = match.group(1), match.group(2)
        out_name = f"{path.stem}.tsv.gz"
        with open(path, "rb") as fin, gzip.open(out_dir / out_name, "wb") as fout:
            shutil.copyfileobj(fin, fout)
        rows.append(
            {
                "contrast_id": f"{target}_vs_{reference}",
                "dataset_key": {"ovary": "elife"}.get(dataset, dataset),
                "variable": "condition",
                "reference": reference,
                "target": target,
                "result_file": out_name,
                "method": "DESeq2",
            }
        )
    manifest = pd.DataFrame(rows)
    manifest.to_csv(out_dir / "contrasts.tsv", sep="\t", index=False)
    print(f"{dataset}: {len(rows)} contrasts -> {out_dir.name}/")


def main() -> None:
    for stale in ("elife_deseq2", "midgut_deseq2"):
        stale_dir = EXPRESSION / stale
        if stale_dir.exists():
            shutil.rmtree(stale_dir)
            print(f"removed superseded {stale}/")
    for dataset in DATASETS:
        samples = bundle_tpm(dataset)
        conditions = sorted({condition_from_sample(s) for s in samples})
        print(f"  conditions: {conditions}")
        bundle_de(dataset)
    shutil.copy(STAGING / "METHODS.md", EXPRESSION / "METHODS.md")
    print("copied METHODS.md")


if __name__ == "__main__":
    main()
