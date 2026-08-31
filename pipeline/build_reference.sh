#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.env
source "$SCRIPT_DIR/config.env"

REFERENCE_ROOT="$RNA_ROOT/reference/$REFERENCE_NAME"
DOWNLOAD_DIR="$REFERENCE_ROOT/downloads"
BUILD_DIR="$REFERENCE_ROOT/build"
INDEX_DIR="$REFERENCE_ROOT/$SALMON_INDEX_NAME"
STATE_DIR="$RNA_ROOT/state"
REFERENCE_COMPLETE="$STATE_DIR/$REFERENCE_COMPLETE_MARKER"
mkdir -p "$DOWNLOAD_DIR" "$BUILD_DIR" "$STATE_DIR"

download_verified() {
  local url="$1"
  local expected_md5="$2"
  local destination="$3"
  local partial="$destination.part"
  local actual_size
  local actual_md5

  if [[ -f "$destination" ]]; then
    actual_md5="$(md5sum "$destination" | awk '{print $1}')"
    if [[ "$actual_md5" != "$expected_md5" ]]; then
      echo "Existing file has unexpected MD5; preserving it and stopping: $destination" >&2
      return 1
    fi
    echo "Already verified: $destination"
    return 0
  fi

  curl \
    --fail \
    --location \
    --retry 12 \
    --retry-delay 5 \
    --retry-all-errors \
    --continue-at - \
    --output "$partial" \
    "$url"

  actual_size="$(wc -c < "$partial")"
  actual_md5="$(md5sum "$partial" | awk '{print $1}')"
  printf 'Downloaded %s bytes; md5=%s\n' "$actual_size" "$actual_md5"
  if [[ "$actual_md5" != "$expected_md5" ]]; then
    echo "MD5 mismatch; preserving partial file and stopping: $partial" >&2
    return 1
  fi
  mv "$partial" "$destination"
}

extract_verified() {
  local archive="$1"
  local expected_md5="$2"
  local destination="$3"
  local partial="$destination.part"
  local actual_md5

  if [[ -f "$destination" ]]; then
    actual_md5="$(md5sum "$destination" | awk '{print $1}')"
    if [[ "$actual_md5" != "$expected_md5" ]]; then
      echo "Existing extracted file has unexpected MD5; preserving it and stopping: $destination" >&2
      return 1
    fi
    echo "Already verified: $destination"
    return 0
  fi
  if [[ -e "$partial" ]]; then
    echo "Partial extraction already exists; preserving it and stopping: $partial" >&2
    return 1
  fi

  pigz --decompress --stdout "$archive" > "$partial"
  actual_md5="$(md5sum "$partial" | awk '{print $1}')"
  if [[ "$actual_md5" != "$expected_md5" ]]; then
    echo "Extracted MD5 mismatch; preserving partial file and stopping: $partial" >&2
    return 1
  fi
  mv "$partial" "$destination"
}

GENOME_GZ="$DOWNLOAD_DIR/$GENOME_FILE"
GTF_GZ="$DOWNLOAD_DIR/$GTF_FILE"
GENOME_FA="$BUILD_DIR/${GENOME_FILE%.gz}"
GENES_GTF="$BUILD_DIR/${GTF_FILE%.gz}"
TRANSCRIPT_RECORDS_GTF="$BUILD_DIR/transcript_records.gtf"
TRANSCRIPTS_FA="$BUILD_DIR/transcripts.filtered.fa"
DECOYS_TXT="$BUILD_DIR/decoys.txt"
GENTROME_FA="$BUILD_DIR/gentrome.fa"
TX2GENE_TSV="$BUILD_DIR/tx2gene.tsv"

download_verified "$REFERENCE_BASE_URL/$GENOME_FILE" "$GENOME_MD5" "$GENOME_GZ"
download_verified "$REFERENCE_BASE_URL/$GTF_FILE" "$GTF_MD5" "$GTF_GZ"
extract_verified "$GENOME_GZ" "$GENOME_UNCOMPRESSED_MD5" "$GENOME_FA"
extract_verified "$GTF_GZ" "$GTF_UNCOMPRESSED_MD5" "$GENES_GTF"

if [[ ! -f "$TRANSCRIPT_RECORDS_GTF" ]]; then
  if [[ -e "$TRANSCRIPT_RECORDS_GTF.part" ]]; then
    echo "Partial filtered GTF exists; preserving it and stopping: $TRANSCRIPT_RECORDS_GTF.part" >&2
    exit 1
  fi
  awk -F '\t' 'BEGIN {OFS = FS} /^#/ || ($3 != "gene" && $9 !~ /transcript_id ""/) {print}' \
    "$GENES_GTF" > "$TRANSCRIPT_RECORDS_GTF.part"
  test -s "$TRANSCRIPT_RECORDS_GTF.part"
  mv "$TRANSCRIPT_RECORDS_GTF.part" "$TRANSCRIPT_RECORDS_GTF"
fi

if [[ ! -f "$TRANSCRIPTS_FA" ]]; then
  if [[ -e "$TRANSCRIPTS_FA.part" ]]; then
    echo "Partial transcriptome exists; preserving it and stopping: $TRANSCRIPTS_FA.part" >&2
    exit 1
  fi
  gffread "$TRANSCRIPT_RECORDS_GTF" -g "$GENOME_FA" -w "$TRANSCRIPTS_FA.part"
  test -s "$TRANSCRIPTS_FA.part"
  mv "$TRANSCRIPTS_FA.part" "$TRANSCRIPTS_FA"
fi

if [[ ! -f "$DECOYS_TXT" ]]; then
  grep '^>' "$GENOME_FA" | cut -d ' ' -f 1 | sed 's/^>//' > "$DECOYS_TXT.part"
  test -s "$DECOYS_TXT.part"
  mv "$DECOYS_TXT.part" "$DECOYS_TXT"
fi

if [[ ! -f "$GENTROME_FA" ]]; then
  if [[ -e "$GENTROME_FA.part" ]]; then
    echo "Partial gentrome exists; preserving it and stopping: $GENTROME_FA.part" >&2
    exit 1
  fi
  cat "$TRANSCRIPTS_FA" "$GENOME_FA" > "$GENTROME_FA.part"
  test -s "$GENTROME_FA.part"
  mv "$GENTROME_FA.part" "$GENTROME_FA"
fi

if [[ ! -f "$TX2GENE_TSV" ]]; then
  python3 "$SCRIPT_DIR/make_tx2gene.py" "$GENES_GTF" "$TX2GENE_TSV"
fi

if [[ ! -f "$REFERENCE_COMPLETE" ]]; then
  if [[ -e "$INDEX_DIR" ]]; then
    echo "Salmon index directory exists without a completion marker; preserving it and stopping: $INDEX_DIR" >&2
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
    printf 'genome_md5=%s\n' "$(md5sum "$GENOME_FA" | awk '{print $1}')"
    printf 'gtf_md5=%s\n' "$(md5sum "$GENES_GTF" | awk '{print $1}')"
    printf 'transcripts=%s\n' "$(grep -c '^>' "$TRANSCRIPTS_FA")"
    printf 'decoys=%s\n' "$(wc -l < "$DECOYS_TXT")"
    du -sh "$REFERENCE_ROOT"
    printf 'keep_duplicates=yes\n'
  } > "$REFERENCE_COMPLETE"
fi

cat "$REFERENCE_COMPLETE"
