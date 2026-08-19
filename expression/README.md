# Gene-expression matrices

These are gene-by-sample TPM (transcripts per million) matrices from our ovary and midgut Salmon analyses plus published neurotranscriptome supplements. TPM is normalized expression abundance: a larger value means that gene contributed a larger share of the sample's sequenced transcript pool. These values are not raw integer read counts.

## Files

The atlas separates **paper** matrices (published values) from **reprocessed**
matrices (raw reads run by us through one identical pipeline; see `METHODS.md`).

### Reprocessed by us — STAR + Salmon, nf-core/rnaseq 3.14.0

- `ovary_star_salmon_gene_tpm.tsv.gz`: 19,920 genes x 33 ovary samples, reprocessed from the Venkataraman et al. (2023) raw FASTQs (`PRJNA796320`).
- `midgut_star_salmon_gene_tpm.tsv.gz`: 19,920 genes x 24 midgut samples from the Vosshall lab midgut raw reads.
- `crop_star_salmon_gene_tpm.tsv.gz`: 19,920 genes x 3 crop samples from the Vosshall lab crop raw reads. A single non-blood-fed condition, so no contrasts exist.
- `ovary_deseq2/`: all 55 pairwise DESeq2 result tables across the eleven ovary conditions.
- `midgut_deseq2/`: all 28 pairwise DESeq2 result tables across the eight midgut conditions.
- `METHODS.md`: the pipeline, reference, parameters, and software versions shared by every reprocessed dataset. The app renders this as its Methods page; it is the single source of truth for reprocessing methodology.

Every reprocessed matrix has `gene_id` and `gene_name` as its first two columns;
every remaining column is one sample's TPM. Sample condition, sex, and replicate
are parsed from the sample names (`Fe.Ov.12hBF.1_S10`), not from a metadata file.

### Published paper values

- `elife_80489_tpm.tsv.gz`: the 18,473-gene published ovary TPM supplement from Venkataraman et al. (2023). Displayed as "Ovary (paper)".
- `neurotranscriptome_2016_aaegl_ru_tpm.tsv.gz`: 16,176 genes x 122 tissue samples using the paper's updated `AaegL.RU` annotation. The first three columns are identifiers/display names.
- `neurotranscriptome_2016_aaegl_3_3_tpm.tsv.gz`: 17,478 genes x 122 tissue samples using the older `AaegL3.3` annotation. The first column is the gene identifier.
- `elife_80489_samples.tsv`: reproductive state, replicate, BioSample accession, and GEO alias for the 33 eLife samples.
- `neurotranscriptome_2016_samples.tsv`: condition, sex, tissue, read length, and mapping metadata from the older paper's library-statistics supplement.
- `neurotranscriptome_2016_gene_annotations.tsv`: paper gene families, OrthoDB categories, Drosophila orthologs/BLASTX hits, and naming evidence for the AaegL.RU genes.

Each differential directory contains `contrasts.tsv` plus the compressed result
tables. The result columns are DESeq2's `baseMean`, `log2FoldChange`, `lfcSE`,
`pvalue`, and `padj`; no differential statistics are computed in the app.

Rebuild the reprocessed bundles from an rsynced results tree with:

```bash
../scripts/bundle_reprocessed_results.py
```

For a new UI, use the `AaegL.RU` matrix as the primary representation of the 2016 paper because it is the authors' updated annotation. Retain `AaegL3.3` as an alternate identifier system for compatibility with older mosquito resources.

Rebuild the published supplementary extracts with:

```bash
../scripts/extract_tpm_matrices.py
```

Rebuild the displayed ovary matrix from a validated Salmon output tree with:

```bash
../scripts/build_ovary_salmon_matrix.py \
  /path/to/PRJNA796320_salmon_results \
  elife_80489_samples.tsv \
  elife_80489_salmon_gene_tpm.tsv.gz
```
