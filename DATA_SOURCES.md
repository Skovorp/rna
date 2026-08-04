# Data sources

The repository includes processed TPM expression matrices and metadata from two open-access studies:

1. Venkataraman K et al. (2023). *Two novel, tightly linked, and rapidly evolving genes underlie Aedes aegypti mosquito reproductive resilience during drought.* eLife 12:e80489. DOI: [10.7554/eLife.80489](https://doi.org/10.7554/eLife.80489). GEO `GSE193470`; BioProject `PRJNA796320`; Zenodo `7758401`. Article and Zenodo materials are CC BY 4.0.
2. Matthews BJ et al. (2016). *The neurotranscriptome of the Aedes aegypti mosquito.* BMC Genomics 17:32. DOI: [10.1186/s12864-015-2239-0](https://doi.org/10.1186/s12864-015-2239-0). BioProject `PRJNA236239`. Article supplementary materials are CC BY 4.0.

It also includes the Nadav Shai / Vosshall lab midgut RNA-seq dataset: 24 paired-end biological libraries spanning non-blood-fed male midgut and female midgut at non-blood-fed, 3, 6, 12, 24, 48, and 72 hours post-blood-meal. The bundled TPM matrix was generated with `nf-core/rnaseq` 3.26.0 and Salmon against the AaegL5 VectorBase 58 + Jové et al. 2019 annotation.

The files under `expression/` contain faithful tabular extracts of the published supplementary tables and the validated gene-level Salmon matrix for the midgut dataset. TPM values are descriptive normalized abundance, not raw read counts.
