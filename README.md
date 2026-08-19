# Aedes RNA Atlas

Small Streamlit prototype for exploring *Aedes aegypti* RNA-seq TPM matrices.

## What it does

- Search individual genes and aliases such as `Ir25a`, `Orco`, `AaegOr7`, and `AAEL005776`, then sort matched genes by each study's mean TPM and toggle them in or out of every result.
- Open the matching genes directly inside the UCSC *Aedes aegypti* Mosquito Cell Atlas, with a gene-colored UMAP followed by a multi-gene expression dot plot. The UMAP can optionally add a second, equally sized plot underneath colored by sex or biological sample. Both embeds have independent atlas-view controls, the dot plot can group nuclei by any categorical atlas annotation, and the atlas authors' default genes can be included or hidden.
- Compare expression across tissues, conditions, and ovary reproductive states, for both published (paper) and reprocessed matrices.
- Browse every pairwise DESeq2 contrast for our reprocessed datasets: 55 ovary contrasts across eleven conditions and 28 midgut contrasts across eight. A configurable FDR threshold colors significant genes in the MA plot. Paper datasets and the single-condition crop dataset are marked `NOT AVAILABLE` rather than being tested from TPM values in the app.
- Read the **Methods** page for the exact pipeline, reference, and parameters shared by every reprocessed dataset, and the **Ovary paper vs reprocessed** page for how closely the published and reprocessed ovary matrices agree.
- Explore IR, OR, GR, and OBP gene families with replicate-aware plots and heatmaps, using one sortable table with a mean-expression column per study to control visibility for predefined and custom families.
- Map biological samples with PCA, UMAP, or t-SNE using all expression genes by default or a smaller most-variable subset.
- Inspect available paper annotations, orthologs, aliases, and raw per-sample TPM values.

## Datasets

| Dataset | Source | Differential expression |
| --- | --- | --- |
| Ovary (paper) | Venkataraman et al., eLife 2023 published TPM | — |
| Ovary (reprocessed) | Our STAR + Salmon run over the same 33 raw samples | 55 pairwise contrasts |
| Atlas (paper) | Matthews et al., BMC Genomics 2016 neurotranscriptome, `AaegL.RU` + legacy `AaegL3.3` | — |
| Midgut (reprocessed) | Vosshall lab midgut RNA-seq | 28 pairwise contrasts |
| Crop (reprocessed) | Vosshall lab crop RNA-seq | — (single condition) |

Reprocessing of the Matthews et al. atlas dataset is still outstanding.

## Run locally

```bash
cd app
./setup.sh
./run.sh
```

Open `http://127.0.0.1:8501`.

## Repository layout

- `app/` — Streamlit app, data layer, and tests.
- `expression/` — our ovary/midgut Salmon TPM matrices, published neurotranscriptome matrices, and compact metadata used by the app.
- `expression/ucsc_mosquito_cell_atlas_genes.json.gz` — compact snapshot of the exact gene names, categorical metadata fields, and curated default genes exposed by all 24 UCSC Mosquito Cell Atlas views; the app uses it only to construct valid embeds.
- `deploy/` — Pi systemd service and one-minute pull/update timer.
- `docs/` — prototype redirect for `rna.getferal.ai`.
- `DATA_SOURCES.md` — paper and dataset attribution.

## Pi prototype

`aedes-rna-atlas-update.timer` checks `main` every minute, runs tests after a fast-forward pull, and restarts the service only when tests pass.

Public HTTPS origin: `https://pi-rus.tailc1209.ts.net`.

Rockefeller-network direct origin: `http://129.85.166.55:8501`. The direct origin
uses the Pi's current DHCP-assigned Ethernet address and may change if its lease changes.

`rna.getferal.ai` is served by GitHub Pages as a redirect to that stable Pi URL. The GoDaddy DNS zone needs a single `rna` CNAME pointing to `skovorp.github.io`.

Initial Pi setup after cloning the repository:

```bash
./deploy/bootstrap-pi.sh
tailscale funnel --bg --yes 8501
```
