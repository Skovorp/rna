# Salmon FASTQ pipeline

This pipeline quantifies ENA FASTQs against one shared current *Aedes aegypti*
reference so datasets can be compared without mixing annotation versions.

## Reference and tools

- Genome: NCBI RefSeq `GCF_002204515.2`, assembly `AaegL5.0`.
- Annotation: NCBI *Aedes aegypti* Annotation Release 101 GTF. Release 101 is
  the latest NCBI annotation release listed for taxon 7159 as of 2026-08-03.
- Salmon: 2.4.1.
- gffread: 0.12.9.

`build_reference.sh` verifies both compressed and expanded NCBI files by MD5,
extracts transcripts with gffread, and builds a decoy-aware Salmon index from
the transcriptome plus genome. It also writes `tx2gene.tsv` for later gene-level
aggregation. Salmon receives that mapping during quantification and writes both
`quant.sf` (transcripts) and `quant.genes.sf` (genes). Index intermediates and
the normalized reference FASTA are retained. The production index uses
`--keepDuplicates`; otherwise Salmon 2.4.1 collapses exact-sequence duplicate
transcripts and can remove otherwise valid genes from `quant.genes.sf`.

NCBI's GTF includes top-level `gene` records whose required-but-empty
`transcript_id ""` field is rejected by gffread 0.12.9. The reference builder
preserves the original verified GTF and writes a second `transcript_records.gtf`
containing the comments and all non-gene records with valid transcript IDs.
Only that derived file is passed to gffread.

## Input contracts

The downloader accepts a tab-separated file with this header:

```text
project run_accession sample_accession secondary_sample_accession sample_alias filename bytes md5 url
```

The quantifier accepts the existing `*_samples.tsv` format in `manifests/`.
One Salmon result is produced per BioSample, combining all runs belonging to
that sample. Both paired- and single-end libraries are supported. Single-end
libraries use a recorded fragment-length assumption of 300 ± 50 bp.

## Persistent worker layout

```text
/rna/
  code_snapshots/  immutable copies of the code and input manifests used
  logs/            tmux-window logs
  quant/           one Salmon output directory per BioSample
  raw/             preserved FASTQ files and any resumable `.part` files
  reference/       downloaded reference, expanded files, transcriptome, index
  software/        downloaded and unpacked tool releases
  state/           completion markers, worklists, success/failure ledgers
```

`monitor_worker.sh` appends disk, download, reference, and quantification
progress to `/rna/logs/status.tsv` once per minute. It is intended to remain in
its own `tmux` window for the lifetime of a worker.

The scripts never delete an existing file. Downloads use `.part` files and are
resumable. If an existing final file fails validation, or an incomplete output
directory is found, the original is preserved and the process stops or writes
to a timestamped attempt directory.

## Adding another dataset

1. Export the FASTQ file manifest and sample manifest using the schemas above.
2. Copy both into a new timestamped `code_snapshots/` directory on the worker.
3. Run `run_project.sh` with the new project accession and those two manifests.
4. Preserve the files under `/rna/state/`; they are the machine-readable run
   ledger and make retries idempotent.

This first version intentionally uses one worker and sample-by-sample Salmon
quantification. A later scheduler can shard samples across workers without
changing the manifests, output contract, or reference provenance.
