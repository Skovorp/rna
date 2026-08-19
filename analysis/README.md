# eLife ovary TPM comparison

`compare_elife_tpm.py` compares the authors' published TPM matrix with the
independent Salmon reanalysis of the same 33 biological samples. It produces a
local HTML report containing the log-expression error distribution, separate
and joint sample-level PCAs, and a cross-source sample-correlation matrix.

Most identifiers in the paper are NCBI `LOC...` genes, while the reanalysis
uses VectorBase `AAEL...` identifiers. `build_gtf_gene_crosswalk.py` resolves
these conservatively from the NCBI Annotation Release 101 and VectorBase 58 +
Jové GTFs: contigs must have a unique matching length and gene coordinates and
strand must match exactly. Ambiguous and non-exact annotation matches are
excluded.

Run the comparison after building or obtaining the crosswalk:

```bash
.venv/bin/python analysis/compare_elife_tpm.py \
  --crosswalk analysis/results/elife_tpm_comparison/annotation_crosswalk.tsv.gz
```

The report uses `log2(TPM + 1)` for agreement statistics. By default, PCA uses
all one-to-one matched genes with no expression or variability cutoff, then
log-transforms and standardizes each gene across samples. Pass a positive
`--top-variable-genes` value to run a most-variable-gene sensitivity analysis.
