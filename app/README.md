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

- **Home:** browse the study catalog with raw-read source links and complete reprocessed TPM downloads. Methodology follows the table. Its final row opens the password screen; private study rows and downloads appear only after unlocking.
- **Genes:** search symbols, `AAEL...` IDs, internal IDs, or aliases; compare panels such as `Ir25a, Orco`; sort matched genes by a separate mean-TPM column for each selected study and toggle them on or off; inspect a UCSC gene-colored UMAP alone or with an equally sized sex/sample-colored plot underneath, plus a multi-gene single-cell expression dot plot; inspect replicate points, group medians, paper annotations, and raw values.
- **Families:** filter to annotated IR, OR, GR, or OBP genes, or enter a custom gene set with the Genes-page token editor; sort matched family members by a separate mean-TPM column for each selected study and toggle any member on or off; optionally use within-gene z-scores to emphasize relative patterns; export the included family matrix. This is not a family-level statistical test.
- **Differential expression:** choose target and reference conditions independently, then browse the precomputed pairwise DESeq2 result in that direction: 55 ovary, 378 tissue-atlas, 28 midgut, and 66 private fat-body / Malpighian-tubule contrasts. Reversing the selector direction reverses log₂ fold change without recomputing DESeq2; base mean, standard error, raw p-value, and adjusted p-value remain the pipeline values. Studies without bundled pipeline outputs display `NOT AVAILABLE` and are never tested from TPM values in the app.
- **Clusters:** select one study and map biological samples with PCA, UMAP, or t-SNE using standardized values from all log-transformed TPM genes by default, or choose a smaller most-variable subset; color points by available sample metadata.

The primary workflows are available from the persistent menu at the top of every page. The Streamlit sidebar and developer toolbar are hidden from the interface.

Morita (2025) leg and Jové (2020) mouthpart TPM tables are available throughout the TPM explorers. Catalog links open the matching study directly. Basrur (2020) reused Matthews reads for Figure 1G–H and deposited six additional Aedes brain libraries at PRJNA612100. Those new libraries have a separate catalog row marked **Brain (reprocessed, pending)**, without an explorer or download until TPM reprocessing is complete. Older exon-view links still open the existing reprocessed Matthews gene-TPM dataset in Genes.

Gene links use only Goldman Table S1.4 and explicit identifiers supplied in the
bundled papers. This connects `Orco`, `AaegOr7`, and `AAEL005776`, and connects
`ppk317`, `ppk00873`, `AAEL000873`, and `LOC5567199`. No prefix stripping,
coordinate, sequence, or shared-symbol inference establishes identity.

Gene details show VectorBase and NCBI descriptions separately. **Methods**
documents the sources and offers a ZIP containing every accessible gene row,
unmatched rows, original IDs, all searchable aliases, both descriptions,
source-table evidence, and the exact comparison pairs. Ambiguous published
aliases remain visible; private study rows require an unlocked session.

This UI reads normalized gene-by-sample matrices. Additional matrices can be added to `expression_explorer/data.py` using the same dataset contract.
