# Paper-vs-reprocessed TPM comparisons

Both comparisons use **identical original identifiers or identifier links
explicitly supplied in paper tables**. `expression/published_gene_aliases.json.gz`
contains the source assertions, citations, table row numbers and checksums.
It is built from Goldman Table S1.4 and the bundled paper tables with:

```bash
python scripts/import_gene_aliases.py --goldman /path/to/GoldmanVosshallShai_TableS1.xlsx
python scripts/rebuild_ovary_comparison.py
python scripts/rebuild_atlas_comparison.py
```

The rebuild scripts update the comparison data and the app's prepared assets.
The old coordinate-based mapping generator and crosswalk have been removed.
No sequence, coordinate, orthology, or shared-symbol inference supplies pairs.
One-to-many and many-to-one candidates are excluded; every retained pair records
its evidence in `matched_genes.tsv.gz`. Download these files together with the
full app mapping from the Methods page.

The ovary comparison uses all 33 biological samples and the smaller subset of
genes with published identifier support. Its current statistics are not directly
comparable to the former genome-wide coordinate-matched analysis.

The tissue-atlas comparison uses all 122 samples present in both sources. It
continues to sum TPM for the 22 repeated historical IDs in the published matrix
before matching gene rows. The three reprocessed libraries without a published
counterpart remain excluded.

PCA settings and sample selections are unchanged: log2(TPM + 1), the 500 most
variable genes by default, centering without unit-variance scaling, and joint
selection by average within-source variance. Pass `--top-variable-genes 0` to
the underlying analysis script for an all-variable-genes sensitivity analysis.
