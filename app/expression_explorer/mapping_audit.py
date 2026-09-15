"""Export the exact gene rows, aliases, and evidence used by the explorer."""
from __future__ import annotations

import gzip
from io import BytesIO
from itertools import combinations
import json
from pathlib import Path
import re
from zipfile import ZIP_DEFLATED, ZipFile

import pandas as pd

from .catalog import PRIVATE_DATASET_KEYS
from .data import load_datasets
from .gene_aliases import DESCRIPTION_COLUMNS, PublishedAliases


def mapping_source_summary(registry: PublishedAliases, expression_dir: Path) -> pd.DataFrame:
    """Count literal source assertions, without expanding connected gene groups.

    Pair counts ignore case and repeated source rows, as the registry does.
    A pair supported by several papers counts once for each supporting source.
    """
    name_pairs = {
        source: {(identifier.casefold(), name.casefold())
                 for row in registry.records if row["source"] == source
                 for identifier in row["identifiers"] for name in row["names"]
                 if identifier.casefold() != name.casefold()}
        for source in registry.sources
    }
    summary = []
    for source in registry.sources:
        records = [row for row in registry.records if row["source"] == source]
        identifiers = {value.casefold() for row in records for value in row["identifiers"]}
        names = {value.casefold() for row in records for value in row["names"]}
        id_pairs = {tuple(sorted(pair)) for row in records
                    for pair in combinations({value.casefold() for value in row["identifiers"]}, 2)}
        loc_aael = {pair for pair in id_pairs
                    if re.fullmatch(r"aael\d+", pair[0]) and re.fullmatch(r"loc\d+", pair[1])}
        other_sources = set().union(*(pairs for key, pairs in name_pairs.items() if key != source))
        local_loc_links = 0
        if source == "jove_2020":
            # These original matrix links stay local; gene123 is not a global ID.
            original = pd.read_csv(expression_dir / registry.sources[source]["file"],
                                   sep="\t", usecols=["gene_id", "gene_name"], keep_default_na=False)
            local = original[original.gene_id.str.fullmatch(r"gene\d+") &
                             original.gene_name.str.fullmatch(r"LOC\d+")]
            local_loc_links = len(local.drop_duplicates())
        summary.append(dict(
            source=source, source_rows=len(records), identifiers=len(identifiers),
            primary_identifiers=len({row["primary"].casefold() for row in records}),
            identifier_pairs=len(id_pairs), loc_aael_pairs=len(loc_aael),
            other_identifier_pairs=len(id_pairs - loc_aael),
            id_name_pairs=len(name_pairs[source]), distinct_names=len(names),
            id_name_pairs_only_in_this_source=len(name_pairs[source] - other_sources),
            id_name_pairs_also_in_goldman=(len(name_pairs[source] & name_pairs.get("goldman_2025_s1.4", set()))
                                          if source != "goldman_2025_s1.4" else 0),
            ppk_names=sum(name.startswith("ppk") for name in names),
            dataset_local_gene_loc_pairs=local_loc_links,
        ))
    return pd.DataFrame(summary)


def mapping_table(datasets: dict, unlocked: bool = False) -> pd.DataFrame:
    frames = []
    for key, dataset in datasets.items():
        if key in PRIVATE_DATASET_KEYS and not unlocked:
            continue
        columns = ["row_id", "stable_id", "internal_id", "raw_symbol", "display_name",
                   "mapping_identity", "mapping_status", "search_text", "search_normalized",
                   "ambiguous_aliases", *DESCRIPTION_COLUMNS, "mapping_sources", "mapping_evidence"]
        frame = dataset.genes.reindex(columns=columns, fill_value="").copy()
        frame.insert(0, "dataset", key)
        frame.insert(1, "study", dataset.label)
        frame.insert(2, "expression_source", dataset.paper)
        frame.insert(3, "annotation_version", dataset.annotation_version)
        frame.insert(4, "matrix_file", dataset.source_file)
        original_columns = [column for column in dataset.genes.columns if column in {
            "Vectorbase Identifier", "Internal gene ID", "Display name", "gene",
            "IDs", "Symbols", "gene_id", "gene_name"}]
        frame["original_annotation"] = [json.dumps(row, ensure_ascii=False)
                                        for row in dataset.genes[original_columns].to_dict("records")]
        if "gene_name" in dataset.genes:
            frame["raw_symbol"] = dataset.genes["gene_name"]
        # Header occupies row 1 in each bundled TSV; row_id itself is zero based.
        frame.insert(6, "matrix_tsv_row", frame["row_id"] + 2)
        frames.append(frame.rename(columns={"stable_id": "original_gene_id",
                                           "raw_symbol": "original_gene_name",
                                           "search_text": "searchable_aliases",
                                           "search_normalized": "normalized_search_aliases"}))
    return pd.concat(frames, ignore_index=True)


def mapping_archive(expression_dir: Path, unlocked: bool = False) -> bytes:
    datasets = load_datasets(expression_dir)
    registry = PublishedAliases.load(expression_dir)
    output = BytesIO()
    with ZipFile(output, "w", compression=ZIP_DEFLATED) as archive:
        archive.writestr("gene_mapping.tsv", mapping_table(datasets, unlocked).to_csv(sep="\t", index=False))
        source_rows = pd.DataFrame(registry.records)
        for column in ("identifiers", "names"):
            source_rows[column] = source_rows[column].map(lambda value: json.dumps(value, ensure_ascii=False))
        archive.writestr("published_source_rows.tsv", source_rows.to_csv(sep="\t", index=False))
        archive.writestr("mapping_source_summary.tsv",
                         mapping_source_summary(registry, expression_dir).to_csv(sep="\t", index=False))
        archive.writestr("sources.json", json.dumps(registry.sources, indent=2))
        archive.writestr("processing_methods.md", (expression_dir / "PROCESSING_DETAILS.md").read_text())
        for study, folder in (("ovary", "elife_tpm_comparison"), ("neurotranscriptome", "atlas_tpm_comparison")):
            path = expression_dir.parent / "analysis" / "results" / folder / "matched_genes.tsv.gz"
            if path.exists():
                archive.writestr(f"{study}_comparison_pairs.tsv", gzip.decompress(path.read_bytes()))
        archive.writestr("README.txt", """Gene mapping audit

gene_mapping.tsv contains one row per expression-matrix gene, including genes
without a published mapping. Original IDs and expression rows are preserved.
mapping_identity groups identifiers explicitly connected by paper rows; an
empty value means only original identifiers are available. It is not permission
to sum or merge expression rows. Multiple rows can share published identifiers.
searchable_aliases and normalized_search_aliases are the exact search indices.
Internal gene123 IDs are scoped to their own dataset.

published_source_rows.tsv records each assertion and its source table row.
mapping_source_summary.tsv counts literal assertions by source, deduplicated
within each paper (case-insensitive). identifier_pairs connect two ID fields;
id_name_pairs connect an ID field and a name field. These categories can overlap
when a named atlas ID is also supplied as a symbol. Counts are not unique genes
and do not include extra aliases reached through connected identifier groups.
id_name_pairs_only_in_this_source excludes pairs repeated in any other source.
ppk_names counts distinct names beginning with ppk, not unique genes.
Jove's dataset_local_gene_loc_pairs are original gene123/LOC pairs retained only
inside that dataset; they are not assertions in the global source-row bundle.
sources.json supplies citations, source filenames, and SHA-256 checksums.
processing_methods.md contains the full RNA-seq commands and parameters.
The VectorBase and NCBI description columns preserve Goldman S1.4 verbatim.
Blank descriptions mean no description was supplied for that gene.

The comparison-pair files list the exact pairs used in the displayed analyses.
Only identical original IDs or explicit published identifier links are used;
ambiguous one-to-many/many-to-one pairs are excluded. Genomic coordinates,
sequence similarity, orthology and shared symbol text do not establish links.
Search can return several genes sharing an alias; those are kept separate.
Private dataset rows are included only after the session is unlocked.
""")
    return output.getvalue()
