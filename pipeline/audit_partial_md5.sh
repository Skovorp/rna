#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 PROJECT FASTQ_MANIFEST [OUTPUT_TSV]" >&2
  exit 2
fi

PROJECT="$1"
MANIFEST="$2"
RNA_ROOT="${RNA_ROOT:-/rna}"
PARALLELISM="${MD5_AUDIT_JOBS:-4}"
RAW_DIR="$RNA_ROOT/raw/$PROJECT"
OUTPUT_TSV="${3:-$RNA_ROOT/state/$PROJECT.partial_md5_audit.$(date -u +%Y%m%dT%H%M%SZ).tsv}"

[[ -f "$MANIFEST" ]] || {
  echo "Missing manifest: $MANIFEST" >&2
  exit 1
}
[[ -d "$RAW_DIR" ]] || {
  echo "Missing raw directory: $RAW_DIR" >&2
  exit 1
}

mkdir -p "$(dirname "$OUTPUT_TSV")"
export PROJECT RAW_DIR

{
  printf 'checked_at\tproject\tfilename\tbytes\texpected_md5\tactual_md5\tstatus\n'
  while IFS=$'\t' read -r filename expected_md5; do
    [[ -f "$RAW_DIR/$filename.part" ]] || continue
    printf '%s\0%s\0' "$filename" "$expected_md5"
  done < <(
    awk -F '\t' -v project="$PROJECT" \
      'NR > 1 && $1 == project {print $6 "\t" $8}' "$MANIFEST"
  ) | xargs -0 -r -P "$PARALLELISM" -n 2 bash -c '
    filename="$1"
    expected_md5="$2"
    path="$RAW_DIR/$filename.part"
    bytes="$(stat -c %s "$path")"
    actual_md5="$(md5sum -- "$path" | awk "{print \$1}")"
    status="mismatch"
    [[ "$actual_md5" == "$expected_md5" ]] && status="match"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$(date --iso-8601=seconds)" "$PROJECT" "$filename" "$bytes" \
      "$expected_md5" "$actual_md5" "$status"
  ' _
} > "$OUTPUT_TSV"

audited="$(awk 'NR > 1 {n++} END {print n + 0}' "$OUTPUT_TSV")"
matches="$(awk -F '\t' 'NR > 1 && $7 == "match" {n++} END {print n + 0}' "$OUTPUT_TSV")"
mismatches="$(awk -F '\t' 'NR > 1 && $7 == "mismatch" {n++} END {print n + 0}' "$OUTPUT_TSV")"

printf 'output=%s\n' "$OUTPUT_TSV"
printf 'audited=%s\n' "$audited"
printf 'matches=%s\n' "$matches"
printf 'mismatches=%s\n' "$mismatches"

