# Aedes RNA Atlas

Small Streamlit prototype for exploring published *Aedes aegypti* RNA-seq TPM matrices.

## What it does

- Search individual genes and aliases such as `Ir25a`, `Orco`, `AaegOr7`, and `AAEL005776`.
- Compare expression across tissues, conditions, and ovary reproductive states.
- Compare every gene between two conditions with a positive-TPM MA plot using readable base-10 axes for average TPM and the exact A/B fold ratio. A configurable FDR threshold colors significant genes gold above a gray background, with TPM summaries, Welch p-values, and Benjamini–Hochberg FDR.
- Explore IR, OR, GR, and OBP gene families with replicate-aware plots and heatmaps.
- Map biological samples with PCA, UMAP, or t-SNE using the most-variable expression genes.
- Inspect available paper annotations, orthologs, aliases, and raw per-sample TPM values.

## Run locally

```bash
cd app
./setup.sh
./run.sh
```

Open `http://127.0.0.1:8501`.

## Repository layout

- `app/` — Streamlit app, data layer, and tests.
- `expression/` — published TPM matrices and compact metadata used by the app.
- `deploy/` — Pi systemd service and one-minute pull/update timer.
- `docs/` — prototype redirect for `rna.getferal.ai`.
- `DATA_SOURCES.md` — paper and dataset attribution.

## Pi prototype

The Pi runs the public app with `RNA_ATLAS_PUBLIC=1`, which removes the local nf-core upload control. `aedes-rna-atlas-update.timer` checks `main` every minute, runs tests after a fast-forward pull, and restarts the service only when tests pass.

Public origin: `https://pi-rus.tailc1209.ts.net`.

`rna.getferal.ai` is served by GitHub Pages as a redirect to that stable Pi URL. The GoDaddy DNS zone needs a single `rna` CNAME pointing to `skovorp.github.io`.

Initial Pi setup after cloning the repository:

```bash
./deploy/bootstrap-pi.sh
tailscale funnel --bg --yes 8501
```
