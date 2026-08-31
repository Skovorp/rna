#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.env
source "$SCRIPT_DIR/config.env"

REFERENCE_NAME="pi_Nadav_Shai_Midgut_RNAseq"
REFERENCE_ROOT="$RNA_ROOT/reference/$REFERENCE_NAME"
BUILD_DIR="$REFERENCE_ROOT/build"
INDEX_DIR="$REFERENCE_ROOT/salmon_index_keep_duplicates"
STATE_DIR="$RNA_ROOT/state"
COMPLETE_MARKER="$STATE_DIR/$REFERENCE_NAME.reference.complete.txt"

GENOME="$REFERENCE_ROOT/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta"
GTF="$REFERENCE_ROOT/AaegLVP_VB58-Jove19_MT_noS1_geneNames.sorted.gtf"
GENOME_SHA256="fd96bcf4d05fa54b0dc4edeebaf077ab6206c35bda1faa150a392dfdcd0545ec"
GTF_SHA256="0bf20c3fae7f8788e44b56fdbc1e81f5dd502da8c4350d5e4c8757a048a74580"

TRANSCRIPTS_FA="$BUILD_DIR/transcripts.fa"
DECOYS_TXT="$BUILD_DIR/decoys.txt"
GENTROME_FA="$BUILD_DIR/gentrome.fa"
TX2GENE_TSV="$BUILD_DIR/tx2gene.tsv"
GENOME_SEQNAMES="$BUILD_DIR/genome.seqnames.txt"
GTF_SEQNAMES="$BUILD_DIR/gtf.seqnames.txt"

mkdir -p "$BUILD_DIR" "$STATE_DIR"

verify_sha256() {
  local file="$1"
  local expected="$2"
  local actual
  [[ -f "$file" ]] || { echo "Missing source file: $file" >&2; return 1; }
  actual="$(sha256sum "$file" | awk '{print $1}')"
  [[ "$actual" == "$expected" ]] || {
    echo "SHA-256 mismatch; preserving file and stopping: $file" >&2
    return 1
  }
  printf 'Verified SHA-256: %s\n' "$file"
}

write_once() {
  local destination="$1"
  shift
  local partial="$destination.part"
  if [[ -f "$destination" ]]; then
    echo "Already present: $destination"
    return 0
  fi
  if [[ -e "$partial" ]]; then
    echo "Partial output exists; preserving it and stopping: $partial" >&2
    return 1
  fi
  "$@" > "$partial"
  [[ -s "$partial" ]]
  mv "$partial" "$destination"
}

verify_sha256 "$GENOME" "$GENOME_SHA256"
verify_sha256 "$GTF" "$GTF_SHA256"

write_once "$GENOME_SEQNAMES" bash -c \
  'grep "^>" "$1" | cut -d " " -f 1 | sed "s/^>//" | sort -u' _ "$GENOME"
write_once "$GTF_SEQNAMES" bash -c \
  'awk -F "\t" "!/^#/ && NF >= 9 {print \$1}" "$1" | sort -u' _ "$GTF"

missing_seqnames="$(comm -23 "$GTF_SEQNAMES" "$GENOME_SEQNAMES" | head -n 20)"
if [[ -n "$missing_seqnames" ]]; then
  echo "GTF sequence names missing from genome; preserving inputs and stopping:" >&2
  printf '%s\n' "$missing_seqnames" >&2
  exit 1
fi

if [[ ! -f "$TRANSCRIPTS_FA" ]]; then
  [[ ! -e "$TRANSCRIPTS_FA.part" ]] || {
    echo "Partial transcriptome exists; preserving it and stopping: $TRANSCRIPTS_FA.part" >&2
    exit 1
  }
  gffread "$GTF" -g "$GENOME" -w "$TRANSCRIPTS_FA.part"
  [[ -s "$TRANSCRIPTS_FA.part" ]]
  mv "$TRANSCRIPTS_FA.part" "$TRANSCRIPTS_FA"
fi

write_once "$DECOYS_TXT" bash -c \
  'grep "^>" "$1" | cut -d " " -f 1 | sed "s/^>//"' _ "$GENOME"
write_once "$GENTROME_FA" bash -c 'cat "$1" "$2"' _ "$TRANSCRIPTS_FA" "$GENOME"

if [[ ! -f "$TX2GENE_TSV" ]]; then
  python3 "$SCRIPT_DIR/make_tx2gene.py" "$GTF" "$TX2GENE_TSV"
fi

if [[ ! -f "$COMPLETE_MARKER" ]]; then
  if [[ -e "$INDEX_DIR" ]]; then
    echo "Index exists without completion marker; preserving it and stopping: $INDEX_DIR" >&2
    exit 1
  fi
  salmon index \
    --transcripts "$GENTROME_FA" \
    --decoys "$DECOYS_TXT" \
    --index "$INDEX_DIR" \
    --threads "$REFERENCE_THREADS" \
    --keepIntermediate \
    --keepFixedFasta \
    --keepDuplicates

  {
    date --iso-8601=seconds
    salmon --version
    printf 'reference_name=%s\n' "$REFERENCE_NAME"
    printf 'genome_sha256=%s\n' "$GENOME_SHA256"
    printf 'gtf_sha256=%s\n' "$GTF_SHA256"
    printf 'transcripts=%s\n' "$(grep -c '^>' "$TRANSCRIPTS_FA")"
    printf 'genes=%s\n' "$(awk -F '\t' 'NR > 1 {seen[$2] = 1} END {print length(seen)}' "$TX2GENE_TSV")"
    printf 'decoys=%s\n' "$(wc -l < "$DECOYS_TXT")"
    printf 'keep_duplicates=yes\n'
    du -sh "$REFERENCE_ROOT"
  } > "$COMPLETE_MARKER"
fi

cat "$COMPLETE_MARKER"
