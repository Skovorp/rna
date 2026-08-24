# Aedes RNA Atlas

Local Streamlit UI for exploring published and reprocessed mosquito RNA-seq expression matrices.

## Start

```bash
cd app
./setup.sh
./run.sh
```

Open `http://localhost:8501`.

## Main workflows

- **Home:** understand the atlas, the available studies, and the four analysis workflows before opening an explorer.
- **Genes:** search symbols, `AAEL...` IDs, internal IDs, or aliases; compare panels such as `Ir25a, Orco`; sort matched genes by a separate mean-TPM column for each selected study and toggle them on or off; inspect a UCSC gene-colored UMAP alone or with an equally sized sex/sample-colored plot underneath, plus a multi-gene single-cell expression dot plot; inspect replicate points, group medians, paper annotations, and raw values.
- **Families:** filter to annotated IR, OR, GR, or OBP genes, or enter a custom gene set with the Genes-page token editor; sort matched family members by a separate mean-TPM column for each selected study and toggle any member on or off; optionally use within-gene z-scores to emphasize relative patterns; export the included family matrix. This is not a family-level statistical test.
- **Differential expression:** choose target and reference conditions independently, then browse the precomputed pairwise DESeq2 result in that direction: 55 ovary, 378 tissue-atlas, 28 midgut, and 66 private fat-body / Malpighian-tubule contrasts. Reversing the selector direction reverses log₂ fold change without recomputing DESeq2; base mean, standard error, raw p-value, and adjusted p-value remain the pipeline values. Studies without bundled pipeline outputs display `NOT AVAILABLE` and are never tested from TPM values in the app.
- **Clusters:** select one study and map biological samples with PCA, UMAP, or t-SNE using standardized values from all log-transformed TPM genes by default, or choose a smaller most-variable subset; color points by available sample metadata.

The primary workflows are available from the persistent menu at the top of every page. The Streamlit sidebar and developer toolbar are hidden from the interface.

The alias layer maps `Orco`, `AaegOr7`, and `AAEL005776` to the same gene. It also strips the historical `Aaeg` prefix for cross-paper matching such as `AaegIr25a` → `Ir25a`.

This UI reads normalized gene-by-sample matrices. Additional matrices can be added to `expression_explorer/data.py` using the same dataset contract.
