# RNA-seq processing & differential expression — methodology

All datasets on this page (**ovary**, **crop**, **midgut**) were processed with the *identical* pipeline, reference, and parameters described below. The only per-dataset difference is the input samplesheet.

## 1. Quantification — nf-core/rnaseq 3.14.0

Standard [nf-core/rnaseq](https://nf-co.re/rnaseq/3.14.0) workflow (STAR genomic alignment + Salmon quantification), run with Nextflow 23.10.1 and Apptainer containers on Slurm:

```bash
nextflow run nf-core/rnaseq \
  -r 3.14.0 \
  -profile apptainer \
  -c slurm.config \
  -params-file params.json \
  -resume
```

`params.json` (identical for every dataset except `input`/`outdir`):

```json
{
  "input":  "<dataset>/samplesheet.csv",
  "outdir": "<dataset>/results",
  "fasta":  "reference/aedes_aaegl5/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta",
  "gtf":    "reference/aedes_aaegl5/AaegLVP_VB58-Jove19_MT_noS1_geneNames.sorted.gtf",
  "save_reference": true,
  "max_cpus": 8,
  "max_memory": "128.GB",
  "max_time": "72.h"
}
```

- **Reference:** *Aedes aegypti* LVP_AGWG AaegL5 genome (VectorBase release 68) with a VB58/Jove-derived GTF (MT included, gene names patched).
- **Samplesheet:** `sample,fastq_1,fastq_2,strandedness` with `strandedness: auto` (inferred per-sample by the pipeline).
- All defaults otherwise: TrimGalore adapter/quality trimming, STAR alignment, Salmon quantification of the STAR transcriptome BAM, Picard MarkDuplicates, RSeQC/Qualimap/dupRadar QC, MultiQC report.
- Per-dataset QC: `results/multiqc/star_salmon/multiqc_report.html`.
- Gene-level count matrices used downstream: `results/star_salmon/salmon.merged.gene_counts_length_scaled.tsv`.

## 2. Differential expression — DESeq2, all pairwise contrasts

DE was run identically for every dataset with ≥2 conditions (ovary, midgut; crop contains a single condition and therefore has no contrasts), using DESeq2 inside the nf-core biocontainer:

```bash
Rscript pairwise_de.R \
  <dataset>/results/star_salmon/salmon.merged.gene_counts_length_scaled.tsv \
  <dataset>/results/de_pairwise
```

Exact procedure (see `pairwise_de.R` for the source):

1. **Input:** Salmon length-scaled gene counts, rounded to integers.
2. **Condition assignment:** parsed from sample names (`Fe.<Tissue>.<condition>.<replicate>`); every condition has ≥2 biological replicates.
3. **Model:** `DESeqDataSetFromMatrix(..., design = ~condition)`; genes with total count < 10 across all samples removed; standard `DESeq()` (median-of-ratios size factors, dispersion shrinkage, Wald test).
4. **Contrasts:** every pairwise combination of conditions. Log2 fold changes shrunk with `lfcShrink(type = "ashr")` (falling back to unshrunk `results()` if ashr is unavailable in the container).
5. **Multiple testing:** Benjamini–Hochberg adjusted p-values (`padj`); significance threshold used in summaries: `padj < 0.05`.
6. **Outputs per dataset** (`results/de_pairwise/`):
   - `DE_<A>_vs_<B>.tsv` — one file per contrast (`A` vs reference `B`): `gene_id`, `gene_name`, `baseMean`, `log2FoldChange`, `lfcSE`, `pvalue`, `padj`, sorted by `padj`.
   - `pca.pdf` / `pca_data.tsv` — PCA of variance-stabilized (VST, blind) counts.
   - `dds.rds` — the fitted DESeq2 object for reproducibility.

## 3. Software versions

| Component | Version |
| --- | --- |
| nf-core/rnaseq | 3.14.0 |
| Nextflow | 23.10.1 |
| Containers | Apptainer, nf-core biocontainer images |
| Quantification | STAR + Salmon (nf-core 3.14.0 defaults) |
| DE | DESeq2 (nf-core deseq2 biocontainer, mulled-v2-8849acf3…) |
| Genome | AaegL5 / LVP_AGWG, VectorBase release 68 |
