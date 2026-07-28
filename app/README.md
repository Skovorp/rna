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

- **Home:** understand the atlas, the available studies, and the four analysis workflows before opening an explorer.
- **Genes:** search symbols, `AAEL...` IDs, internal IDs, or aliases; compare panels such as `Ir25a, Orco`; inspect replicate points, group medians, paper annotations, and raw values.
- **Families:** filter to annotated IR, OR, GR, or OBP genes, or enter a custom gene set with the Genes-page token editor; rank individual genes by mean TPM across all samples; show all genes by default or select the top N; optionally use within-gene z-scores to emphasize relative patterns; export the complete family matrix. This is not a family-level statistical test.
- **Compare conditions:** compare every gene between two conditions in one study; view an MA plot with readable base-10 axes for average TPM and the exact A/B fold ratio. A configurable FDR threshold colors significant genes gold and draws them above gray nonsignificant genes. The table retains TPM summaries, raw Welch p-values, and Benjamini–Hochberg FDR. Genes with zero mean TPM in either condition remain in the table but are omitted from the plot because their ratio is undefined.
- **Clusters:** select one study and map biological samples with PCA, UMAP, or t-SNE using standardized values from the most-variable log-transformed TPM genes; color points by available sample metadata.

The primary workflows are available from the persistent menu at the top of every page. The Streamlit sidebar and developer toolbar are hidden from the interface.

The alias layer maps `Orco`, `AaegOr7`, and `AAEL005776` to the same gene. It also strips the historical `Aaeg` prefix for cross-paper matching such as `AaegIr25a` → `Ir25a`.

This UI reads normalized gene-by-sample matrices. Additional matrices can be added to `expression_explorer/data.py` using the same dataset contract.
