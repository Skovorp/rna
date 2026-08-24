# Aedes RNA Atlas

Small Streamlit prototype for exploring *Aedes aegypti* RNA-seq TPM matrices.

## What it does

- Search individual genes and aliases such as `Ir25a`, `Orco`, `AaegOr7`, and `AAEL005776`, then sort matched genes by each study's mean TPM and toggle them in or out of every result.
- Open the matching genes directly inside the UCSC *Aedes aegypti* Mosquito Cell Atlas, with a gene-colored UMAP followed by a multi-gene expression dot plot. The UMAP can optionally add a second, equally sized plot underneath colored by sex or biological sample. Both embeds have independent atlas-view controls, the dot plot can group nuclei by any categorical atlas annotation, and the atlas authors' default genes can be included or hidden.
- Compare expression across tissues, conditions, and ovary reproductive states, for both published (paper) and reprocessed matrices.
- Browse every pairwise DESeq2 contrast for our reprocessed datasets: 55 ovary, 378 tissue-atlas, 28 midgut, and 66 private fat-body / Malpighian-tubule contrasts. Separate target and reference selectors expose either fold-change direction without recomputing the underlying DESeq2 test. Paper datasets and the single-condition crop dataset are marked `NOT AVAILABLE` rather than being tested from TPM values in the app.
- Read the **Methods** page for the exact shared pipeline and the **Ovary paper vs reprocessed** and **Tissue atlas paper vs reprocessed** pages for direct comparisons with the published matrices.
- Explore IR, OR, GR, and OBP gene families with replicate-aware plots and heatmaps, using one sortable table with a mean-expression column per study to control visibility for predefined and custom families.
- Map biological samples with PCA, UMAP, or t-SNE using all expression genes by default or a smaller most-variable subset.
- Inspect available paper annotations, orthologs, aliases, and raw per-sample TPM values.

## Datasets

| Dataset | Source | Differential expression |
| --- | --- | --- |
| Ovary (paper) | Venkataraman et al., eLife 2023 published TPM | — |
| Ovary (reprocessed) | Our STAR + Salmon run over the same 33 raw samples | 55 pairwise contrasts |
| Atlas (paper) | Matthews et al., BMC Genomics 2016 neurotranscriptome, `AaegL.RU` + legacy `AaegL3.3` | — |
| Atlas (reprocessed) | Our STAR + Salmon run over the same raw reads | 378 pairwise contrasts |
| Midgut (reprocessed) | Vosshall lab midgut RNA-seq | 28 pairwise contrasts |
| Fat body & Malpighian tubules (reprocessed, private) | Vosshall lab blood-meal time course | 66 pairwise contrasts |
| Crop (reprocessed, private) | Vosshall lab crop RNA-seq | — (single condition) |

The tissue-atlas comparison uses the 122 samples present in the paper matrix.
Three additional recovered libraries remain visible in the reprocessed dataset
but are excluded from the comparison because they have no published counterpart.

## Run locally

```bash
cd app
./setup.sh
./run.sh
```

Open `http://127.0.0.1:8501`.

## Repository layout

- `app/` — Streamlit app, data layer, and tests.
- `expression/` — reprocessed Salmon TPM matrices, published matrices, DESeq2 results, and compact metadata used by the app.
- `expression/ucsc_mosquito_cell_atlas_genes.json.gz` — compact snapshot of the exact gene names, categorical metadata fields, and curated default genes exposed by all 24 UCSC Mosquito Cell Atlas views; the app uses it only to construct valid embeds.
- `deploy/` — systemd, updater, and nginx configuration for production.
- `DATA_SOURCES.md` — paper and dataset attribution.

## Production

Production is `https://mosquito.rockefeller.edu` on the Hetzner VPS. The
`rna-atlas-update.timer` checks `main` every minute, fast-forwards the checkout,
warms the derived dataset cache, and restarts the Streamlit service. Nginx is
reloaded only when its tracked configuration changes and validates.
