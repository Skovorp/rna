"""Study provenance and downloads, available without loading expression matrices."""

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class CatalogEntry:
    key: str
    label: str
    description: str
    raw_sources: tuple[tuple[str, str], ...] = ()
    comparison: str | None = None
    explorer: str | None = None
    tpm_file: str | None = None
    private: bool = False
    pending_reprocessing: bool = False


OVARY_RAW = (("PRJNA796320", "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA796320"),)
NEURO_RAW = (("PRJNA236239", "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA236239"),)

CATALOG = (
    CatalogEntry(
        "ovary_paper", "Ovary (published)",
        "Published TPM from [Venkataraman et al. (2023)](https://www.biorxiv.org/content/10.1101/2022.03.01.482582).",
        OVARY_RAW, "/Ovary_paper_vs_reprocessed",
    ),
    CatalogEntry(
        "elife", "Ovary (reprocessed)",
        "Our STAR + Salmon gene TPM for the same 33 samples, with 55 pairwise DESeq2 contrasts.",
        OVARY_RAW, "/Ovary_paper_vs_reprocessed",
        tpm_file="ovary_star_salmon_gene_tpm.tsv.gz",
    ),
    CatalogEntry(
        "neuro_ru", "Neurotranscriptome (published)",
        "Published AaegL.RU and legacy AaegL3.3 TPM from [Matthews et al. (2016)](https://www.biorxiv.org/content/10.1101/026823).",
        NEURO_RAW, "/Atlas_paper_vs_reprocessed",
    ),
    CatalogEntry(
        "atlas", "Neurotranscriptome (reprocessed)",
        "Our STAR + Salmon gene TPM from the same raw reads, with 378 pairwise DESeq2 contrasts. [Basrur et al. (2020)](https://elifesciences.org/articles/63982) reused these reads for the fruitless exon measurements in Figure 1G–H.",
        NEURO_RAW, "/Atlas_paper_vs_reprocessed",
        tpm_file="atlas_star_salmon_gene_tpm.tsv.gz",
    ),
    CatalogEntry(
        "midgut", "Midgut (reprocessed)",
        "Vosshall lab midgut RNA-seq: our STAR + Salmon gene TPM and 28 pairwise DESeq2 contrasts.",
        tpm_file="midgut_star_salmon_gene_tpm.tsv.gz",
    ),
    CatalogEntry(
        "basrur", "Brain (reprocessed, pending)",
        "[Basrur et al. (2020)](https://elifesciences.org/articles/63982): six Aedes aegypti brain samples (three female, three male).",
        (("PRJNA612100", "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA612100"),),
        pending_reprocessing=True,
    ),
    CatalogEntry(
        "morita", "Legs, wild type & Orco mutants (published)",
        "[Morita et al. (2025)](https://www.science.org/doi/full/10.1126/sciadv.adn5758): gene TPM summed from the authors’ Salmon transcript tables for ten female leg samples.",
        (("PRJNA1020561", "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1020561"),),
        explorer="/?page=Genes&study=morita",
    ),
    CatalogEntry(
        "jove", "Mouthparts (published)",
        "[Jové et al. (2020)](https://www.cell.com/neuron/fulltext/S0896-6273(20)30719-4): published stylet and labium TPM, with four replicates per group.",
        (("PRJNA605870", "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA605870"),),
        explorer="/?page=Genes&study=jove",
    ),
    CatalogEntry(
        "yedlin", "Fat body & Malpighian tubules (reprocessed, private)",
        "Our STAR + Salmon gene TPM over the blood-meal time course, with 66 pairwise DESeq2 contrasts.",
        tpm_file="yedlin_star_salmon_gene_tpm.tsv.gz", private=True,
    ),
    CatalogEntry(
        "crop", "Crop (reprocessed, private)",
        "Our STAR + Salmon gene TPM from three non-blood-fed replicates. One condition, so no differential contrasts.",
        tpm_file="crop_star_salmon_gene_tpm.tsv.gz", private=True,
    ),
)

CATALOG_BY_KEY = {entry.key: entry for entry in CATALOG}
PRIVATE_DATASET_KEYS = frozenset(entry.key for entry in CATALOG if entry.private)
DATASET_LABELS = {entry.key: entry.label for entry in CATALOG}
DATASET_LABELS["neuro_ru"] = "Neurotranscriptome (published, AaegL.RU)"
DATASET_LABELS["neuro_legacy"] = "Neurotranscriptome (published, legacy AaegL3.3)"


def visible_catalog(unlocked: bool = False) -> tuple[CatalogEntry, ...]:
    return tuple(entry for entry in CATALOG if unlocked or not entry.private)


def tpm_download_path(key: str, expression_dir: Path, unlocked: bool = False) -> Path:
    """Resolve only declared downloads, checking access before reading a file."""
    entry = CATALOG_BY_KEY[key]
    if entry.private and not unlocked:
        raise PermissionError("Unlock private datasets before downloading their tables.")
    if entry.tpm_file is None:
        raise ValueError(f"No reprocessed TPM download for {key}")
    return expression_dir / entry.tpm_file
