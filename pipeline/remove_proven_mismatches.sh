#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 PROJECT AUDIT_TSV" >&2
  exit 2
fi

PROJECT="$1"
AUDIT_TSV="$2"
RNA_ROOT="${RNA_ROOT:-/rna}"
RAW_DIR="$RNA_ROOT/raw/$PROJECT"
STATE_DIR="$RNA_ROOT/state"

[[ "${ALLOW_DELETE_PROVEN_MISMATCHES:-no}" == "yes" ]] || {
  echo "Refusing deletion without ALLOW_DELETE_PROVEN_MISMATCHES=yes" >&2
  exit 1
}
[[ -f "$AUDIT_TSV" ]] || {
  echo "Missing audit TSV: $AUDIT_TSV" >&2
  exit 1
}
[[ -d "$RAW_DIR" ]] || {
  echo "Missing raw directory: $RAW_DIR" >&2
  exit 1
}

if pgrep -f "^bash .*/download_project[.]sh $PROJECT " >/dev/null || \
   pgrep -f "^curl .*${RAW_DIR}/" >/dev/null; then
  echo "Refusing deletion while the project downloader is active" >&2
  exit 1
fi

timestamp="$(date -u +%Y%m%dT%H%M%SZ)"
ledger="$STATE_DIR/$PROJECT.proven_mismatch_deletion.$timestamp.tsv"
record="$STATE_DIR/$PROJECT.proven_mismatch_deletion.$timestamp.txt"
mkdir -p "$STATE_DIR"

awk -F '\t' '
  NR == 1 {next}
  NF != 7 {bad = 1; next}
  $7 != "mismatch" {bad = 1}
  seen[$3]++ {bad = 1}
  END {exit bad ? 1 : 0}
' "$AUDIT_TSV" || {
  echo "Audit contains malformed, duplicate, or non-mismatch rows" >&2
  exit 1
}

rows="$(awk 'NR > 1 {n++} END {print n + 0}' "$AUDIT_TSV")"
(( rows > 0 )) || {
  echo "Audit contains no mismatch rows" >&2
  exit 1
}

printf 'deleted_at\tproject\tpath\tbytes\texpected_md5\taudited_actual_md5\treason\n' > "$ledger"
total_bytes=0
while IFS=$'\t' read -r checked_at project filename bytes expected_md5 actual_md5 status; do
  [[ "$project" == "$PROJECT" && "$status" == "mismatch" ]] || exit 1
  [[ "$expected_md5" =~ ^[0-9a-f]{32}$ && "$actual_md5" =~ ^[0-9a-f]{32}$ ]] || exit 1
  [[ "$expected_md5" != "$actual_md5" ]] || exit 1
  path="$RAW_DIR/$filename.part"
  [[ -f "$path" ]] || {
    echo "Audited partial is missing: $path" >&2
    exit 1
  }
  current_bytes="$(stat -c %s "$path")"
  [[ "$current_bytes" == "$bytes" ]] || {
    echo "Audited partial size changed: $path" >&2
    exit 1
  }
  total_bytes=$((total_bytes + current_bytes))
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$(date --iso-8601=seconds)" "$PROJECT" "$path" "$current_bytes" \
    "$expected_md5" "$actual_md5" "user_authorized_proven_md5_mismatch" >> "$ledger"
done < <(tail -n +2 "$AUDIT_TSV")

{
  printf 'prepared_at=%s\n' "$(date --iso-8601=seconds)"
  printf 'project=%s\n' "$PROJECT"
  printf 'authorization=delete_provably_broken_files\n'
  printf 'audit=%s\n' "$AUDIT_TSV"
  printf 'audit_sha256=%s\n' "$(sha256sum "$AUDIT_TSV" | awk '{print $1}')"
  printf 'ledger=%s\n' "$ledger"
  printf 'ledger_sha256_before_deletion=%s\n' "$(sha256sum "$ledger" | awk '{print $1}')"
  printf 'files=%s\n' "$rows"
  printf 'bytes=%s\n' "$total_bytes"
} > "$record"
sync "$ledger" "$record"

while IFS=$'\t' read -r deleted_at project path bytes expected_md5 actual_md5 reason; do
  rm -- "$path"
done < <(tail -n +2 "$ledger")

remaining=0
while IFS=$'\t' read -r deleted_at project path bytes expected_md5 actual_md5 reason; do
  [[ ! -e "$path" ]] || remaining=$((remaining + 1))
done < <(tail -n +2 "$ledger")
(( remaining == 0 )) || {
  echo "$remaining authorized mismatch files still exist" >&2
  exit 1
}

{
  printf 'deleted_at=%s\n' "$(date --iso-8601=seconds)"
  printf 'deleted_files=%s\n' "$rows"
  printf 'deleted_bytes=%s\n' "$total_bytes"
  printf 'ledger_sha256_after_deletion=%s\n' "$(sha256sum "$ledger" | awk '{print $1}')"
} >> "$record"
sync "$record"

cat "$record"
