# Aedes RNA Atlas

Local Streamlit UI for exploring the published TPM matrices from the two mosquito papers in this workspace.

## Start

```bash
cd app
./setup.sh
./run.sh
```

Open `http://localhost:8501`.

## Main workflows

- **Genes:** search symbols, `AAEL...` IDs, internal IDs, or aliases; compare panels such as `Ir25a, Orco`; inspect replicate points, group medians, paper annotations, and raw values.
- **Families:** inspect IR, OR, GR, or OBP families across separate study panels; rank by peak group median; view relative-pattern heatmaps; export the complete family matrix.
- **Data & provenance:** use the collapsed sidebar for dataset descriptions, interpretation limits, and local nf-core/rnaseq imports.

The updated AaegL.RU annotation is the default. The AaegL3.3 matrix remains available for legacy identifier compatibility but is not presented as a separate biological experiment.

The import dialog accepts `salmon.merged.gene_tpm.tsv`, `rsem.merged.gene_tpm.tsv`, or equivalent nf-core/rnaseq 3.26.0 merged gene-TPM output. Imported matrices remain local under `rna/expression/imports/`.

The alias layer maps `Orco`, `AaegOr7`, and `AAEL005776` to the same gene. It also strips the historical `Aaeg` prefix for cross-paper matching such as `AaegIr25a` → `Ir25a`.

nf-core/rnaseq remains the upstream route for future raw-read processing. This UI reads normalized gene-by-sample matrices; additional nf-core-derived matrices can be added to `expression_explorer/data.py` using the same dataset contract.
