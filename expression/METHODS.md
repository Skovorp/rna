## RNA-seq processing

**Quantification:** all completed reprocessed datasets use nf-core/rnaseq 3.14.0 (Nextflow 23.10.1, Apptainer), TrimGalore, STAR + Salmon, and automatic strandedness. Reference: AaegL5 / VectorBase 68 genome with the VB58/Jové GTF, mitochondrial genes included and gene names patched. Published datasets retain the authors' expression values.

**Differential expression:** DESeq2 on rounded Salmon length-scaled counts; genes with total count <10 removed; design `~condition`; all pairwise contrasts with ≥2 biological replicates per condition. Wald tests, Benjamini–Hochberg `padj <0.05`, and ashr fold-change shrinkage when available. PCA uses blind variance-stabilized counts. Full commands and parameters are in the mapping ZIP.

## Gene matching

We use identical original IDs and links explicitly provided in Goldman **Table S1.4** or the paper tables below. ID links work both ways: **NCBI `LOC…` ↔ VectorBase `AAEL…`**. Names are searchable aliases; a shared name alone does not join genes. Comparisons include only unambiguous one-to-one ID pairs. Unlinked genes keep their original IDs.

**VectorBase description** and **NCBI description** are copied separately from Goldman S1.4; missing descriptions stay blank.
