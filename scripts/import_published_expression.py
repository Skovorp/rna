"""Extract published tables without changing their units or sample membership.

Requires pandas and openpyxl for this one-time build; neither Excel nor source
archives are read by the app. Source URLs and checksums are written with outputs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import zipfile

import numpy as np
import pandas as pd


MORITA_COMMIT = "a8a1fdcc11491967fae9cbb593e18aebd930cb63"
JOVE_COMMIT = "78d9b066340dd2f25e8f8b40d4fec7a8f99652f3"
SOURCES = {
    "morita_orco_vs_wt_quants.zip": f"https://raw.githubusercontent.com/VosshallLab/Morita_Vosshall2023/{MORITA_COMMIT}/figure8/orco_vs_wt_quants.zip",
    "morita_annotation.gtf.zip": f"https://raw.githubusercontent.com/VosshallLab/Morita_Vosshall2023/{MORITA_COMMIT}/figure8/AaegLVP_VB58-Jove19.gtf.zip",
    "morita_chemoreceptors.tsv": f"https://raw.githubusercontent.com/VosshallLab/Morita_Vosshall2023/{MORITA_COMMIT}/figure8/chemoreceptor_gene_list.tsv",
    "jove_2020_tpm.csv": f"https://raw.githubusercontent.com/VosshallLab/Jove_Vosshall_2020/{JOVE_COMMIT}/RNAseq_merged_annotation/merge_19_2_3_TPM_final.csv",
    "jove_data_file_1.xlsx": f"https://raw.githubusercontent.com/VosshallLab/Jove_Vosshall_2020/{JOVE_COMMIT}/Data_File_1_Neuron.xlsx",
}


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(
        path, sep="\t", index=False, float_format="%.17g",
        compression={"method": "gzip", "mtime": 0} if path.suffix == ".gz" else None,
    )


def validate_values(frame: pd.DataFrame) -> None:
    values = frame.to_numpy(dtype=float)
    assert np.isfinite(values).all(), "Non-finite published expression values"
    assert (values >= 0).all(), "Negative published expression values"
    assert frame.index.is_unique and frame.columns.is_unique


def build_morita(source: Path, out: Path) -> dict:
    with zipfile.ZipFile(source / "morita_annotation.gtf.zip") as archive:
        gtf = archive.read("AaegLVP_VB58-Jove19.gtf").decode()
    tx_to_gene = {}
    for line in gtf.splitlines():
        columns = line.split("\t")
        if len(columns) != 9:
            continue
        attrs = dict(re.findall(r'(\w+) "([^"]+)"', columns[8]))
        if "transcript_id" in attrs and "gene_id" in attrs:
            tx, gene = attrs["transcript_id"], attrs["gene_id"]
            assert tx not in tx_to_gene or tx_to_gene[tx] == gene
            tx_to_gene[tx] = gene
    gene_series = []
    rows = []
    source_totals = {}
    with zipfile.ZipFile(source / "morita_orco_vs_wt_quants.zip") as archive:
        for group in ("WT", "KO"):
            for replicate in range(1, 6):
                sample = f"{group}{replicate}"
                with archive.open(f"orco_vs_wt_quants/{sample}_quant/quant.sf") as handle:
                    quant = pd.read_csv(handle, sep="\t").set_index("Name")
                validate_values(quant[["TPM"]])
                genes = quant.index.to_series().map(tx_to_gene)
                assert genes.notna().all(), f"Unmapped transcripts in {sample}"
                # tximport's gene abundance is the sum of transcript TPM for
                # each gene. Preserve that abundance; do not normalize counts.
                summed = quant["TPM"].groupby(genes).sum().rename(sample)
                assert np.isclose(summed.sum(), quant["TPM"].sum(), rtol=1e-12)
                gene_series.append(summed)
                source_totals[sample] = float(quant["TPM"].sum())
                condition = "Wild type" if group == "WT" else "Orco mutant (5/16)"
                rows.append(dict(sample=sample, sex="female", tissue="forelegs and midlegs",
                                 genotype=condition, condition=group, condition_label=condition,
                                 tissue_condition=f"Legs / {condition}", replicate=replicate,
                                 reproductive_state="Not specified"))
    values = pd.concat(gene_series, axis=1)
    validate_values(values)
    symbols = pd.read_csv(source / "morita_chemoreceptors.tsv", sep="\t").set_index("id")["name"]
    assert symbols.index.is_unique
    frame = values.copy()
    frame.insert(0, "gene_name", frame.index.to_series().map(symbols).fillna(""))
    frame.insert(0, "gene_id", frame.index)
    write_tsv(frame, out / "morita_2025_gene_tpm.tsv.gz")
    write_tsv(pd.DataFrame(rows), out / "morita_2025_samples.tsv")
    return {"genes": len(values), "samples": len(rows), "unit": "TPM",
            "aggregation": "Sum author Salmon transcript TPM by the authors' GTF gene_id (tximport abundance).",
            "transcript_tpm_sums": source_totals}


def build_jove(source: Path, out: Path) -> dict:
    original = pd.read_csv(source / "jove_2020_tpm.csv").set_index("GeneID")
    annotated = pd.read_excel(source / "jove_data_file_1.xlsx", sheet_name="TPM_all_coding_transcripts").set_index("GeneID")
    name_column = annotated.columns[0]
    assert set(original.index) == set(annotated.index)
    assert original.index.is_unique and annotated.index.is_unique
    values = original.copy()
    validate_values(values)
    workbook_values = annotated.loc[values.index, values.columns].to_numpy(dtype=float)
    # The workbook rounds the higher-precision CSV (largest discrepancy <5e-7
    # TPM). Keep the CSV values, and verify every zero remains a true zero.
    np.testing.assert_allclose(values, workbook_values, rtol=0, atol=5e-7)
    np.testing.assert_array_equal(values.to_numpy() == 0, workbook_values == 0)
    names = annotated.loc[values.index, name_column].astype(str)
    # Preserve the paper's gene and LOC identifiers. Named receptors retain
    # their exact published symbols and are resolved against atlas aliases.
    frame = values.copy()
    frame.insert(0, "gene_name", names)
    frame.insert(0, "gene_id", values.index)
    write_tsv(frame, out / "jove_2020_gene_tpm.tsv.gz")
    rows = []
    for sample in values.columns:
        group, replicate = re.fullmatch(r"(Female|Labium|Male)(\d+)", sample).groups()
        sex = "male" if group == "Male" else "female"
        tissue = "labium" if group == "Labium" else "stylet"
        condition = f"{sex.capitalize()} {tissue}"
        rows.append(dict(sample=sample, sex=sex, tissue=tissue, condition=group,
                         condition_label=condition, tissue_condition=condition, replicate=int(replicate),
                         reproductive_state="Not specified"))
    write_tsv(pd.DataFrame(rows), out / "jove_2020_samples.tsv")
    return {"genes": len(values), "samples": len(rows), "unit": "TPM",
            "validation": "All CSV values match Data File 1 within its 5e-7 TPM rounding tolerance; zero patterns match exactly.",
            "max_workbook_rounding_difference": float(np.abs(values.to_numpy() - workbook_values).max())}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("expression"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    report = {"sources": {name: {"url": url, "sha256": hashlib.sha256((args.sources / name).read_bytes()).hexdigest()}
                          for name, url in SOURCES.items()}}
    for key, builder in (("morita", build_morita), ("jove", build_jove)):
        report[key] = builder(args.sources, args.output)
        print(key, report[key], flush=True)
    (args.output / "published_2020_2025_provenance.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
