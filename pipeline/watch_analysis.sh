#!/usr/bin/env bash
set -euo pipefail

RNA_ROOT="${RNA_ROOT:-/rna}"
INTERVAL_SECONDS="${INTERVAL_SECONDS:-1800}"
LOG="$RNA_ROOT/logs/analysis_watchdog.tsv"
REFERENCE_MARKER="$RNA_ROOT/state/pi_Nadav_Shai_Midgut_RNAseq.reference.complete.txt"
QUANT_MARKER="$RNA_ROOT/state/PRJNA796320.quant.complete.txt"
QUANT_SUCCESS="$RNA_ROOT/state/PRJNA796320.quant.success.tsv"
QUANT_FAILURE="$RNA_ROOT/state/PRJNA796320.quant.failure.tsv"
ATLAS_MARKER="$RNA_ROOT/state/PRJNA236239.download.complete.txt"

mkdir -p "$RNA_ROOT/logs"
if [[ ! -f "$LOG" ]]; then
  printf 'checked_at\tassessment\treference_complete\treference_index_active\tquant_complete\tquant_success\tquant_failures\tquant_active\tatlas_download_complete\tatlas_downloader_active\tatlas_curls\tatlas_finals\tatlas_partials\tdisk_available_bytes\n' > "$LOG"
fi

while true; do
  reference_complete=0
  reference_index_active="$(pgrep -fc '[s]almon index' || true)"
  quant_complete=0
  quant_success=0
  quant_failures=0
  quant_active="$(pgrep -fc '[s]almon quant' || true)"
  atlas_download_complete=0
  atlas_downloader_active="$(pgrep -fc '[d]ownload_project.sh PRJNA236239' || true)"
  atlas_curls="$(ps -eo args | grep '[c]url .*\/rna\/raw\/PRJNA236239\/' | wc -l)"
  atlas_finals="$(find "$RNA_ROOT/raw/PRJNA236239" -maxdepth 1 -type f -name '*.fastq.gz' | wc -l)"
  atlas_partials="$(find "$RNA_ROOT/raw/PRJNA236239" -maxdepth 1 -type f -name '*.part' | wc -l)"
  disk_available_bytes="$(df -B1 "$RNA_ROOT" | awk 'NR == 2 {print $4}')"

  [[ -f "$REFERENCE_MARKER" ]] && reference_complete=1
  [[ -f "$QUANT_MARKER" ]] && quant_complete=1
  [[ -f "$ATLAS_MARKER" ]] && atlas_download_complete=1
  if [[ -f "$QUANT_SUCCESS" ]]; then
    quant_success="$(awk -F '\t' 'NR > 1 {seen[$3] = 1} END {print length(seen)}' "$QUANT_SUCCESS")"
  fi
  if [[ -f "$QUANT_FAILURE" ]]; then
    # Count only failed samples that do not yet have a later successful run.
    # The failure ledger is append-only, so recovered attempts remain auditable.
    quant_failures="$(
      awk -F '\t' '
        FILENAME == ARGV[1] && FNR > 1 {successful[$3] = 1; next}
        FILENAME == ARGV[2] && FNR > 1 && !($3 in successful) {failed[$3] = 1}
        END {print length(failed)}
      ' "$QUANT_SUCCESS" "$QUANT_FAILURE"
    )"
  fi

  assessment="ok"
  if (( disk_available_bytes < 50 * 1024 * 1024 * 1024 )); then
    assessment="alert_disk_below_50GiB"
  elif (( quant_failures > 0 && quant_active == 0 )); then
    assessment="alert_quant_failure"
  elif (( reference_complete == 0 && reference_index_active == 0 )); then
    assessment="alert_reference_stalled"
  elif (( reference_complete == 1 && quant_complete == 0 && quant_success < 33 && quant_active == 0 )); then
    assessment="alert_quant_stalled"
  elif (( atlas_download_complete == 0 && atlas_downloader_active == 0 )); then
    assessment="alert_atlas_download_stalled"
  fi

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$(date --iso-8601=seconds)" "$assessment" "$reference_complete" \
    "$reference_index_active" "$quant_complete" "$quant_success" \
    "$quant_failures" "$quant_active" "$atlas_download_complete" \
    "$atlas_downloader_active" "$atlas_curls" "$atlas_finals" \
    "$atlas_partials" "$disk_available_bytes" >> "$LOG"

  sleep "$INTERVAL_SECONDS"
done
