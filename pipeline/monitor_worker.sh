#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.env
source "$SCRIPT_DIR/config.env"

INTERVAL_SECONDS="${INTERVAL_SECONDS:-60}"
STATUS_FILE="$RNA_ROOT/logs/status.tsv"
mkdir -p "$RNA_ROOT/logs"

if [[ ! -f "$STATUS_FILE" ]]; then
  printf 'timestamp\tdisk_used_bytes\tdisk_available_bytes\tPRJNA796320_raw_bytes\tPRJNA236239_raw_bytes\tfinal_fastqs\tpartial_fastqs\treference_complete\tPRJNA796320_quants\tPRJNA236239_quants\tactive_curls\tactive_salmon\n' \
    > "$STATUS_FILE"
fi

while true; do
  timestamp="$(date --iso-8601=seconds)"
  disk_used="$(df -PB1 "$RNA_ROOT" | awk 'NR == 2 {print $3}')"
  disk_available="$(df -PB1 "$RNA_ROOT" | awk 'NR == 2 {print $4}')"
  small_raw="$(du -sb "$RNA_ROOT/raw/PRJNA796320" 2>/dev/null | awk '{print $1 + 0}')"
  large_raw="$(du -sb "$RNA_ROOT/raw/PRJNA236239" 2>/dev/null | awk '{print $1 + 0}')"
  final_fastqs="$(find "$RNA_ROOT/raw" -type f -name '*.fastq.gz' 2>/dev/null | wc -l)"
  partial_fastqs="$(find "$RNA_ROOT/raw" -type f -name '*.part' 2>/dev/null | wc -l)"
  reference_complete=0
  [[ -f "$RNA_ROOT/state/$REFERENCE_COMPLETE_MARKER" ]] && reference_complete=1
  small_quants="$(find "$RNA_ROOT/quant/PRJNA796320" -type f -name 'quant.complete.txt' 2>/dev/null | wc -l)"
  large_quants="$(find "$RNA_ROOT/quant/PRJNA236239" -type f -name 'quant.complete.txt' 2>/dev/null | wc -l)"
  active_curls="$(pgrep -fc 'curl.*ftp.sra.ebi' || true)"
  active_salmon="$(pgrep -fc 'salmon (index|quant)' || true)"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$timestamp" "$disk_used" "$disk_available" "$small_raw" "$large_raw" \
    "$final_fastqs" "$partial_fastqs" "$reference_complete" "$small_quants" \
    "$large_quants" "$active_curls" "$active_salmon" | tee -a "$STATUS_FILE"
  sleep "$INTERVAL_SECONDS"
done
