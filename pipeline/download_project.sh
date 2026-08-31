#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 PROJECT FASTQ_FILES_TSV" >&2
  exit 2
fi

PROJECT="$1"
MANIFEST="$2"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.env
source "$SCRIPT_DIR/config.env"

DESTINATION="$RNA_ROOT/raw/$PROJECT"
STATE_DIR="$RNA_ROOT/state"
ATTEMPT="$(date -u +%Y%m%dT%H%M%SZ)"
WORKLIST="$STATE_DIR/$PROJECT.download.$ATTEMPT.tsv"
mkdir -p "$DESTINATION" "$STATE_DIR"

if [[ ! -f "$MANIFEST" ]]; then
  echo "Missing FASTQ manifest: $MANIFEST" >&2
  exit 1
fi

awk -F '\t' -v project="$PROJECT" 'NR > 1 && $1 == project' "$MANIFEST" > "$WORKLIST"
if [[ ! -s "$WORKLIST" ]]; then
  echo "No FASTQs found for $PROJECT in $MANIFEST" >&2
  exit 1
fi

required_bytes="$(awk -F '\t' '{total += $7} END {printf "%.0f", total}' "$WORKLIST")"
present_bytes=0
while IFS=$'\t' read -r project run sample secondary alias filename expected_bytes expected_md5 url; do
  output="$DESTINATION/$filename"
  partial="$output.part"
  if [[ -f "$output" ]]; then
    actual_bytes="$(wc -c < "$output")"
  elif [[ -f "$partial" ]]; then
    actual_bytes="$(wc -c < "$partial")"
  else
    actual_bytes=0
  fi
  if (( actual_bytes > expected_bytes )); then
    actual_bytes="$expected_bytes"
  fi
  present_bytes=$((present_bytes + actual_bytes))
done < "$WORKLIST"
remaining_bytes=$((required_bytes - present_bytes))
available_kib="$(df -Pk "$DESTINATION" | awk 'NR == 2 {print $4}')"
available_bytes=$((available_kib * 1024))
printf 'project=%s files=%s required_bytes=%s present_bytes=%s remaining_bytes=%s available_bytes=%s\n' \
  "$PROJECT" "$(wc -l < "$WORKLIST")" "$required_bytes" "$present_bytes" \
  "$remaining_bytes" "$available_bytes"
if (( available_bytes < remaining_bytes + 20 * 1024 * 1024 * 1024 )); then
  echo "Not enough disk for $PROJECT plus a 20 GiB reserve" >&2
  exit 1
fi

download_one() {
  local row="$1"
  local project run sample secondary alias filename expected_bytes expected_md5 url
  local output partial actual_bytes actual_md5
  IFS=$'\t' read -r project run sample secondary alias filename expected_bytes expected_md5 url <<< "$row"
  output="$DESTINATION/$filename"
  partial="$output.part"

  if [[ -f "$output" ]]; then
    actual_bytes="$(wc -c < "$output")"
    actual_md5="$(md5sum "$output" | awk '{print $1}')"
    if [[ "$actual_bytes" == "$expected_bytes" && "$actual_md5" == "$expected_md5" ]]; then
      echo "Already verified: $filename"
      return 0
    fi
    echo "Existing final file failed validation; preserving it and stopping: $output" >&2
    return 1
  fi

  if [[ -f "$partial" ]]; then
    actual_bytes="$(wc -c < "$partial")"
    if (( actual_bytes > expected_bytes )); then
      echo "Partial file is larger than expected; preserving it and stopping: $partial" >&2
      return 1
    fi
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

  actual_bytes="$(wc -c < "$partial")"
  actual_md5="$(md5sum "$partial" | awk '{print $1}')"
  if [[ "$actual_bytes" != "$expected_bytes" || "$actual_md5" != "$expected_md5" ]]; then
    echo "Integrity check failed; preserving partial file: $partial" >&2
    return 1
  fi
  mv "$partial" "$output"
  echo "Verified: $filename ($run / $sample)"
}

export DESTINATION
export -f download_one
xargs -P "$DOWNLOAD_JOBS" -I '{}' bash -c 'download_one "$1"' _ '{}' < "$WORKLIST"

while IFS=$'\t' read -r project run sample secondary alias filename expected_bytes expected_md5 url; do
  output="$DESTINATION/$filename"
  [[ -f "$output" ]] || { echo "Missing downloaded file: $output" >&2; exit 1; }
  actual_bytes="$(wc -c < "$output")"
  actual_md5="$(md5sum "$output" | awk '{print $1}')"
  [[ "$actual_bytes" == "$expected_bytes" && "$actual_md5" == "$expected_md5" ]] || {
    echo "Final validation failed: $output" >&2
    exit 1
  }
done < "$WORKLIST"

{
  date --iso-8601=seconds
  printf 'project=%s\n' "$PROJECT"
  printf 'manifest=%s\n' "$MANIFEST"
  printf 'manifest_md5=%s\n' "$(md5sum "$MANIFEST" | awk '{print $1}')"
  printf 'files=%s\n' "$(wc -l < "$WORKLIST")"
  printf 'bytes=%s\n' "$required_bytes"
  du -sh "$DESTINATION"
} > "$STATE_DIR/$PROJECT.download.complete.txt"

cat "$STATE_DIR/$PROJECT.download.complete.txt"
