"""Gene links asserted in paper tables; no coordinate or symbol inference."""
from __future__ import annotations

from collections import defaultdict
import gzip
import json
from pathlib import Path
import re


DESCRIPTION_COLUMNS = ("vectorbase_AAEL", "vectorbase_symbol", "vectorbase_description",
                       "ncbi_symbol", "ncbi_description")


def identifier_key(value: str) -> str:
    # Punctuation is significant for identity, even though the search UI is
    # forgiving. In particular, do not invent version/prefix conversions.
    return value.strip().casefold()


class PublishedAliases:
    def __init__(self, payload: dict):
        if payload.get("policy") != "explicit_paper_rows_only":
            raise ValueError("Gene aliases must come from explicit paper rows")
        self.sources = payload["sources"]
        self.records = payload["records"]
        parents = {}

        def root(key):
            parents.setdefault(key, key)
            if parents[key] != key:
                parents[key] = root(parents[key])
            return parents[key]

        # Join only identifiers explicitly listed together, never two records
        # just because their symbols/descriptions resemble or equal each other.
        for record in self.records:
            ids = [identifier_key(v) for v in record["identifiers"]]
            for key in ids[1:]:
                parents[root(key)] = root(ids[0])
            root(ids[0])
        self.identities = {key: root(key) for key in parents}
        self.groups = defaultdict(list)
        for record in self.records:
            self.groups[self.identities[identifier_key(record["identifiers"][0])]].append(record)
        self.alias_groups = defaultdict(set)
        self.group_aliases = {}
        for key, records in self.groups.items():
            aliases = sorted({v for r in records for v in r["identifiers"] + r["names"] if v})
            self.group_aliases[key] = tuple(aliases)
            for alias in aliases:
                self.alias_groups[identifier_key(alias)].add(key)

    @classmethod
    def load(cls, expression_dir: Path | str):
        with gzip.open(Path(expression_dir) / "published_gene_aliases.json.gz", "rt") as handle:
            return cls(json.load(handle))

    def identity(self, identifier: str) -> str | None:
        if re.fullmatch(r"gene\d+", identifier, re.I):
            return None
        return self.identities.get(identifier_key(identifier))

    def describe(self, identifier: str, raw_symbol: str = "") -> dict:
        # Some papers place their NCBI identifier in the name column beside a
        # paper-local gene123 ID. Only an explicit stable ID may bridge it.
        anchor = identifier
        key = self.identity(anchor)
        if key is None and re.fullmatch(r"(?:LOC|AAEL)\d+", raw_symbol):
            anchor = raw_symbol
            key = self.identity(anchor)
        if key is None:
            return dict(aliases=(), display_name=raw_symbol or identifier,
                        mapping_identity="", mapping_status="original identifiers only",
                        mapping_sources="", mapping_evidence="", ambiguous_aliases="",
                        **{column: "" for column in DESCRIPTION_COLUMNS})
        records = self.groups[key]
        goldman = [r for r in records if r["source"] == "goldman_2025_s1.4"]
        preferred = goldman or [r for r in records if r["source"] == "morita_2025"] or records
        direct = [r for r in preferred if identifier_key(r["primary"]) == identifier_key(anchor)]
        names = {r["preferred"] for r in direct or preferred if r["preferred"]}
        display = next(iter(names)) if len(names) == 1 else raw_symbol or identifier
        ambiguous = [a for a in self.group_aliases[key] if len(self.alias_groups[identifier_key(a)]) > 1]
        # Keep different genes sharing a published name separate in the UI.
        if display and len(self.alias_groups[identifier_key(display)]) > 1:
            stable_ids = sorted({v for r in records for v in r["identifiers"] if re.fullmatch(r"AAEL\d+", v)})
            display = f"{display} ({stable_ids[0] if len(stable_ids) == 1 else identifier})"
        return dict(
            aliases=self.group_aliases[key], display_name=display, mapping_identity=key,
            mapping_status="ambiguous published aliases" if ambiguous or len({r['primary'] for r in goldman}) > 1 else "published mapping",
            mapping_sources=" | ".join(sorted({r["source"] for r in records})),
            mapping_evidence=json.dumps([dict(source=r["source"], row=r["row"], identifiers=r["identifiers"], names=r["names"]) for r in records], ensure_ascii=False),
            ambiguous_aliases=" | ".join(ambiguous),
            **{column: " | ".join(dict.fromkeys(r.get(column, "") for r in goldman if r.get(column, ""))) for column in DESCRIPTION_COLUMNS},
        )
