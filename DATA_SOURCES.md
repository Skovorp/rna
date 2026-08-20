# Data sources

The repository includes processed TPM expression matrices and metadata from two open-access studies:

1. Venkataraman K et al. (2023). *Two novel, tightly linked, and rapidly evolving genes underlie Aedes aegypti mosquito reproductive resilience during drought.* eLife 12:e80489. DOI: [10.7554/eLife.80489](https://doi.org/10.7554/eLife.80489). GEO `GSE193470`; BioProject `PRJNA796320`; Zenodo `7758401`. Article and Zenodo materials are CC BY 4.0.
2. Matthews BJ et al. (2016). *The neurotranscriptome of the Aedes aegypti mosquito.* BMC Genomics 17:32. DOI: [10.1186/s12864-015-2239-0](https://doi.org/10.1186/s12864-015-2239-0). BioProject `PRJNA236239`. Article supplementary materials are CC BY 4.0.

It also includes the Nadav Shai / Vosshall lab midgut RNA-seq dataset: 24 paired-end biological libraries spanning non-blood-fed male midgut and female midgut at non-blood-fed, 3, 6, 12, 24, 48, and 72 hours post-blood-meal. The bundled TPM matrix was generated with `nf-core/rnaseq` 3.26.0 and Salmon against the AaegL5 VectorBase 58 + Jové et al. 2019 annotation.

The Genes page also embeds the [UCSC Aedes aegypti Mosquito Cell Atlas](https://cells.ucsc.edu/?ds=mosquito+all), from Goldman OV et al. (2025), *A single-nucleus transcriptomic atlas of the adult Aedes aegypti mosquito*, Cell 188:7267–7290.e26, DOI [10.1016/j.cell.2025.10.008](https://doi.org/10.1016/j.cell.2025.10.008). This is a deep-linked external visualization, not a locally reprocessed expression dataset. Its values are normalized single-nucleus expression rather than TPM.

The files under `expression/` contain the validated gene-level Salmon matrices for the ovary and midgut datasets, faithful tabular extracts of the 2016 study's published supplementary tables, and every pairwise DESeq2 contrast produced by `nf-core/differentialabundance` 2.0.0: 28 contrasts across the eight midgut conditions and 55 across the eleven ovary conditions. TPM values are descriptive normalized abundance, not raw read counts; the app displays differential-expression statistics only from the precomputed count-aware pipeline outputs.

`expression/ucsc_mosquito_cell_atlas_genes.json.gz` is a compact routing manifest derived from each UCSC view's public `exprMatrix.json` and `dataset.json`. It records which exact identifiers can be passed as `gene=` for all 24 leaf datasets, plus the categorical metadata fields and author-curated default genes used to configure the embedded multi-gene dot plot. Refresh it with `scripts/update_ucsc_cell_atlas_manifest.py`; no UCSC expression matrix is downloaded or analyzed by the app.

## Displayed-data provenance

The atlas separates **paper** datasets (published values, shown as the authors
released them) from **reprocessed** datasets (raw reads run through our own
pipeline). Where both exist for the same samples, a comparison page reports how
closely they agree.

| Atlas dataset | Displayed values | Differential expression |
| --- | --- | --- |
| Ovary (paper) | Venkataraman et al. published TPM supplement | Not available |
| Ovary (reprocessed) | Our STAR + Salmon gene TPM from all 33 `PRJNA796320` raw samples | All 55 pairwise DESeq2 contrasts |
| Atlas (paper), AaegL.RU | Matthews et al. published `AaegL.RU` TPM matrix | Not available |
| Atlas (paper), legacy AaegL3.3 | Matthews et al. published legacy matrix, retained for identifier compatibility | Not available |
| Midgut (reprocessed) | Our STAR + Salmon gene TPM from the Vosshall lab midgut raw reads | All 28 pairwise DESeq2 contrasts |
| Crop (reprocessed) | Our STAR + Salmon gene TPM from the Vosshall lab crop raw reads | Not applicable — a single condition, so no contrasts exist |

Reprocessing of the Matthews et al. neurotranscriptome (the "atlas" dataset) is
still outstanding; only the published version is shown.

Every reprocessed dataset above went through the *identical* pipeline,
reference, and parameters. The exact workflow, versions, and parameter files
are recorded in `expression/METHODS.md`, which the app renders as its Methods
page. That file is the single source of truth for reprocessing methodology --
do not restate pipeline parameters elsewhere.

## Paper-vs-reprocessed comparisons

`app/assets/ovary_comparison/` holds the self-contained comparison reports
rendered by the atlas:

- `elife_ovary_tpm_full_report.html` -- per-sample and per-gene agreement
  between the published and reprocessed ovary TPM matrices.
- `elife_ovary_zero_nonzero_transitions.html` -- genes that are exactly zero in
  one matrix but expressed in the other.

The generating analysis, its summary statistics, and the identifier crosswalk
live under `analysis/results/elife_tpm_comparison/`.
