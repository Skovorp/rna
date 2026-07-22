# Gene-expression matrices

These are the authors' published gene-by-sample TPM (transcripts per million) matrices. TPM is normalized expression abundance: a larger value means that gene contributed a larger share of the sample's sequenced transcript pool. These values are not raw integer read counts.

## Files

- `elife_80489_tpm.tsv.gz`: 18,473 genes × 33 ovary samples from Venkataraman et al. (2023). The first two columns are `IDs` and `Symbols`; every remaining column is a sample's TPM value.
- `neurotranscriptome_2016_aaegl_ru_tpm.tsv.gz`: 16,176 genes × 122 tissue samples using the paper's updated `AaegL.RU` annotation. The first three columns are identifiers/display names.
- `neurotranscriptome_2016_aaegl_3_3_tpm.tsv.gz`: 17,478 genes × 122 tissue samples using the older `AaegL3.3` annotation. The first column is the gene identifier.
- `elife_80489_samples.tsv`: reproductive state and replicate metadata for the 33 eLife samples.
- `neurotranscriptome_2016_samples.tsv`: condition, sex, tissue, read length, and mapping metadata from the older paper's library-statistics supplement.
- `neurotranscriptome_2016_gene_annotations.tsv`: paper gene families, OrthoDB categories, Drosophila orthologs/BLASTX hits, and naming evidence for the AaegL.RU genes.

For a new UI, use the `AaegL.RU` matrix as the primary representation of the 2016 paper because it is the authors' updated annotation. Retain `AaegL3.3` as an alternate identifier system for compatibility with older mosquito resources.

Rebuild these files with:

```bash
../scripts/extract_tpm_matrices.py
```
