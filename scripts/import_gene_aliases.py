#!/usr/bin/env python3
"""Bundle literal gene identifiers/names from supplied paper tables."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re

import openpyxl


def clean(value):
    return "" if value is None else str(value).strip()


def build(workbook: Path, expression: Path) -> dict:
    records = []
    sources = {}

    def source(key, path, citation, table):
        sources[key] = dict(file=path.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                            citation=citation, table=table)

    def add(source_key, row, identifiers, names, primary, preferred="", descriptions=None):
        identifiers = list(dict.fromkeys(clean(v) for v in identifiers if clean(v)))
        names = list(dict.fromkeys(clean(v) for v in names if clean(v)))
        if not identifiers:
            return
        records.append(dict(source=source_key, row=row, identifiers=identifiers,
                            names=names, primary=primary, preferred=preferred,
                            **(descriptions or {})))

    sheet_name = "Table S1.4 - Gene Look Up"
    source("goldman_2025_s1.4", workbook, "https://doi.org/10.1016/j.cell.2025.10.008", sheet_name)
    book = openpyxl.load_workbook(workbook, read_only=True, data_only=True)
    iterator = book[sheet_name].iter_rows(values_only=True)
    columns = next(iterator)
    for number, values in enumerate(iterator, 2):
        row = {k: clean(v) for k, v in zip(columns, values) if k}
        primary = row["Aedesaegypti_Mosquito_Cell_Atlas_geneID"]
        if not primary:
            continue
        vectorbase = row["vectorbase_AAEL"]
        identifiers = [primary, row["ncbi_symbol"]]
        if vectorbase != "manual annotation":
            identifiers.append(vectorbase)
        add("goldman_2025_s1.4", number, identifiers,
            [row["Aedesaegypti_Mosquito_Cell_Atlas_genesymbol"], row["vectorbase_symbol"]],
            primary, row["Aedesaegypti_Mosquito_Cell_Atlas_genesymbol"],
            {k: row[k] for k in ("vectorbase_AAEL", "vectorbase_symbol", "vectorbase_description",
                                 "ncbi_symbol", "ncbi_description")})
    book.close()

    configurations = [
        ("matthews_2016_s5", "neurotranscriptome_2016_aaegl_ru_tpm.tsv.gz",
         "Vectorbase Identifier", "Display name", "https://doi.org/10.1186/s12864-015-2239-0", "Additional File 5"),
        ("morita_2025", "morita_2025_gene_tpm.tsv.gz", "gene_id", "gene_name",
         "https://doi.org/10.1126/sciadv.adn5758", "Author chemoreceptor_gene_list.tsv / gene TPM extract"),
        ("venkataraman_2023", "elife_80489_tpm.tsv.gz", "IDs", "Symbols",
         "https://doi.org/10.7554/eLife.80489", "Supplementary Data File 02"),
        ("jove_2020", "jove_2020_gene_tpm.tsv.gz", "gene_id", "gene_name",
         "https://doi.org/10.1016/j.neuron.2020.09.019", "Data File 1 / TPM CSV"),
    ]
    for key, filename, id_column, name_column, citation, table in configurations:
        path = expression / filename
        source(key, path, citation, table)
        with gzip.open(path, "rt") as handle:
            for number, row in enumerate(csv.DictReader(handle, delimiter="\t"), 2):
                identifier, name = clean(row[id_column]), clean(row[name_column])
                # Internal gene123 identifiers belong to one paper only. They
                # remain searchable locally, never become cross-paper links.
                if re.fullmatch(r"gene\d+", identifier) or not name or name == identifier:
                    continue
                add(key, number, [identifier], [name], identifier, name)
    return dict(schema_version=1, policy="explicit_paper_rows_only", sources=sources, records=records)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--goldman", required=True, type=Path)
    parser.add_argument("--expression", type=Path, default=Path(__file__).resolve().parents[1] / "expression")
    args = parser.parse_args()
    payload = build(args.goldman, args.expression)
    output = args.expression / "published_gene_aliases.json.gz"
    output.write_bytes(gzip.compress(json.dumps(payload, ensure_ascii=False, sort_keys=True).encode(), mtime=0))
    print(f"{len(payload['records']):,} source rows -> {output}")


if __name__ == "__main__":
    main()
