"""Bundle complete source TPMs and the frozen mappings used by each comparison."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import gzip
import hashlib
from io import BytesIO
import json
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

import numpy as np
import pandas as pd

from scripts.rebuild_atlas_comparison import normalize_sample
from scripts.rebuild_ovary_comparison import published_sample_key, reprocessed_sample_key


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class ComparisonInputs:
    published: str
    reprocessed: str
    metadata: str
    results: str
    assets: str
    published_id: str
    published_annotation_columns: int
    paper: str
    accession: str


COMPARISONS = {
    "ovary": ComparisonInputs(
        "elife_80489_tpm.tsv.gz", "ovary_star_salmon_gene_tpm.tsv.gz",
        "elife_80489_samples.tsv", "elife_tpm_comparison", "ovary_comparison",
        "IDs", 2, "https://doi.org/10.7554/eLife.80489", "PRJNA796320",
    ),
    "neurotranscriptome": ComparisonInputs(
        "neurotranscriptome_2016_aaegl_ru_tpm.tsv.gz", "atlas_star_salmon_gene_tpm.tsv.gz",
        "neurotranscriptome_2016_samples.tsv", "atlas_tpm_comparison", "atlas_comparison",
        "Vectorbase Identifier", 3, "https://doi.org/10.1186/s12864-015-2239-0", "PRJNA236239",
    ),
}


def sample_matching(study: str, published: list[str], reprocessed: list[str]) -> pd.DataFrame:
    """Use the comparison rebuilders' existing sample-name rules."""
    published_key, reprocessed_key = ((published_sample_key, reprocessed_sample_key)
                                      if study == "ovary" else (normalize_sample, normalize_sample))
    lookup = {reprocessed_key(sample): sample for sample in reprocessed}
    if len(lookup) != len(reprocessed) or len({published_key(s) for s in published}) != len(published):
        raise ValueError("Ambiguous sample names in comparison inputs")
    rows = []
    used = set()
    for order, sample in enumerate(published, 1):
        counterpart = lookup[published_key(sample)]
        used.add(counterpart)
        rows.append(dict(comparison_order=order, published_sample=sample,
                         reprocessed_sample=counterpart, included_in_comparison=True, exclusion_reason=""))
    for sample in reprocessed:
        if sample not in used:
            rows.append(dict(comparison_order="", published_sample="", reprocessed_sample=sample,
                             included_in_comparison=False, exclusion_reason="No published TPM counterpart"))
    return pd.DataFrame(rows)


def build_download(study: str, root: Path = ROOT) -> Path:
    spec = COMPARISONS[study]
    expression = root / "expression"
    results = root / "analysis" / "results" / spec.results
    assets = root / "app" / "assets" / spec.assets
    files = {
        "published_tpm.tsv": gzip.decompress((expression / spec.published).read_bytes()),
        "reprocessed_tpm.tsv": gzip.decompress((expression / spec.reprocessed).read_bytes()),
        "gene_matching.tsv": gzip.decompress((results / "matched_genes.tsv.gz").read_bytes()),
        "sample_metadata.tsv": (expression / spec.metadata).read_bytes(),
        "analysis_summary.json": (results / "summary.json").read_bytes(),
        "mapping_sources.json": json.dumps(json.loads(gzip.decompress(
            (expression / "published_gene_aliases.json.gz").read_bytes()))["sources"], indent=2).encode(),
    }
    published = pd.read_csv(BytesIO(files["published_tpm.tsv"]), sep="\t", keep_default_na=False)
    reprocessed = pd.read_csv(BytesIO(files["reprocessed_tpm.tsv"]), sep="\t", keep_default_na=False)
    genes = pd.read_csv(BytesIO(files["gene_matching.tsv"]), sep="\t", keep_default_na=False)
    published_samples = published.columns.tolist()[spec.published_annotation_columns:]
    reprocessed_samples = reprocessed.columns.tolist()[2:]
    samples = sample_matching(study, published_samples, reprocessed_samples)
    files["sample_matching.tsv"] = samples.to_csv(sep="\t", index=False).encode()
    summary = json.loads(files["analysis_summary.json"])
    # Refuse to publish data for different figures or an outdated gene mapping.
    if summary != json.loads((assets / "figures.json").read_text())["summary"]:
        raise ValueError("Comparison summary and live figures differ; rebuild the comparison first")
    if genes.published_id.duplicated().any() or genes.reanalysis_id.duplicated().any():
        raise ValueError("Gene matching must be one-to-one")
    if reprocessed.gene_id.duplicated().any():
        raise ValueError("Duplicate reprocessed gene IDs")
    if study == "ovary" and published[spec.published_id].duplicated().any():
        raise ValueError("Duplicate published ovary gene IDs")
    included = samples[samples.included_in_comparison]
    # This aggregation is the existing Matthews protocol, not a new gene link.
    collapsed = published.groupby(spec.published_id, sort=False)[published_samples].sum()
    x = collapsed.loc[genes.published_id, included.published_sample].to_numpy(dtype=float)
    y = reprocessed.set_index("gene_id").loc[genes.reanalysis_id, included.reprocessed_sample].to_numpy(dtype=float)
    if (len(genes), len(included), len(collapsed), len(reprocessed)) != (
        summary["matched_genes"], summary["samples"], summary["published_genes"], summary["reanalysis_genes"]
    ):
        raise ValueError("Input dimensions differ from the displayed comparison")
    if not (np.allclose(x.mean(axis=1), genes.published_mean_tpm, rtol=1e-10, atol=1e-10)
            and np.allclose(y.mean(axis=1), genes.reanalysis_mean_tpm, rtol=1e-10, atol=1e-10)):
        raise ValueError("TPM inputs differ from the frozen gene matching; rebuild the comparison first")

    duplicates = len(published) - len(collapsed)
    excluded = len(samples) - len(included)
    manifest = dict(
        schema_version=1, comparison=study, paper=spec.paper, raw_read_accession=spec.accession,
        published_id_column=spec.published_id, reprocessed_id_column="gene_id",
        published_duplicate_handling="sum TPM by original identifier" if duplicates else "no duplicate identifiers",
        published_rows=len(published), published_unique_ids=len(collapsed), reprocessed_rows=len(reprocessed),
        published_samples=len(published_samples), reprocessed_samples=len(reprocessed_samples),
        compared_samples=len(included), excluded_reprocessed_samples=excluded, matched_genes=len(genes),
        source_files={
            filename: dict(path=f"expression/{original}", sha256=hashlib.sha256((expression / original).read_bytes()).hexdigest())
            for filename, original in (("published_tpm.tsv", spec.published),
                                       ("reprocessed_tpm.tsv", spec.reprocessed),
                                       ("sample_metadata.tsv", spec.metadata))
        },
    )
    files["README.md"] = f"""# {study.capitalize()}: comparison source data

These are the complete original TPM tables, including unmatched genes and all
source sample columns. They are byte-for-byte decompressions of the bundled
source files. No gene rows or values have been renamed, filtered, or summed in
these downloads. The tables contain TPM, not sequencing reads or raw counts.

## Contents

- `published_tpm.tsv`: {len(published):,} original rows × {len(published_samples)} samples.
- `reprocessed_tpm.tsv`: {len(reprocessed):,} genes × {len(reprocessed_samples)} samples.
- `gene_matching.tsv`: the {len(genes):,} one-to-one pairs used in the analysis,
  with mapping methods, source-table evidence, per-gene statistics and PCA flags.
- `sample_matching.tsv`: original sample columns in comparison order;
  `included_in_comparison` selects the {len(included)} compared samples.
  {excluded} extra reprocessed samples are explicitly marked as excluded.
- `sample_metadata.tsv`: original published sample annotations. Its `sample`
  (ovary) or `library_id` (neurotranscriptome) matches `published_sample`.
- `analysis_summary.json`: statistics and PCA settings shown on the website.
- `mapping_sources.json`: source filenames, citations and checksums for mapping evidence.
- `manifest.json`: input provenance, dimensions and SHA-256 checksums.

## Align the tables

Use the exported gene and sample matching; do not match names or column position.
The published ID column is `{spec.published_id}`; the reprocessed ID is `gene_id`.
There are {duplicates} extra published rows with repeated identifiers. Sum TPM
for each repeated original ID **before** selecting matched genes, as the analysis
does. Keep the complete original rows when inspecting the source data.

```python
import json
import numpy as np
import pandas as pd

manifest = json.load(open("manifest.json"))
published = pd.read_csv("published_tpm.tsv", sep="\\t")
reprocessed = pd.read_csv("reprocessed_tpm.tsv", sep="\\t")
genes = pd.read_csv("gene_matching.tsv", sep="\\t")
samples = pd.read_csv("sample_matching.tsv", sep="\\t")
samples = samples[samples.included_in_comparison].sort_values("comparison_order")
published = published.groupby(manifest["published_id_column"], sort=False)[samples.published_sample].sum()
x = published.loc[genes.published_id, samples.published_sample].to_numpy(float)
y = reprocessed.set_index("gene_id").loc[genes.reanalysis_id, samples.reprocessed_sample].to_numpy(float)
xp, yp = np.log2(x + 1), np.log2(y + 1)
print("Pearson r:", np.corrcoef(xp.ravel(), yp.ravel())[0, 1])
print("Median absolute log2 error:", np.median(np.abs(yp - xp)))
```

For PCA, use the per-panel `used_for_*_pca` flags in `gene_matching.tsv`,
log2(TPM + 1), and centered genes without unit-variance scaling. The summary
records the full selection protocol. Sample identity uses all matched genes.
Unmatched genes stay in the source tables but are excluded from comparison
statistics. Identical original IDs or explicit published links establish pairs;
ambiguous links are excluded. No new gene mapping is performed by this export.

Paper: {spec.paper}
Reprocessed reads: https://www.ncbi.nlm.nih.gov/bioproject/{spec.accession}
Analysis code: https://github.com/Skovorp/rna/tree/main/analysis
""".encode()
    manifest["files_sha256"] = {name: hashlib.sha256(data).hexdigest() for name, data in files.items()}
    files["manifest.json"] = json.dumps(manifest, indent=2, sort_keys=True).encode()
    destination = assets / f"{study}_comparison_data.zip"
    temporary = destination.with_suffix(".zip.tmp")
    try:
        with ZipFile(temporary, "w", compression=ZIP_DEFLATED, compresslevel=6) as archive:
            for name, data in sorted(files.items()):
                info = ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
                info.compress_type = ZIP_DEFLATED
                archive.writestr(info, data)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    print(f"{study}: {len(genes):,} gene pairs, {len(included)} compared samples; {destination.stat().st_size:,} bytes", flush=True)
    return destination


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", nargs="?", choices=COMPARISONS)
    args = parser.parse_args()
    for name in ([args.study] if args.study else COMPARISONS):
        build_download(name)
