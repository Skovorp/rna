#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 PROJECT FASTQ_FILES_TSV PROJECT_SAMPLES_TSV" >&2
  exit 2
fi

PROJECT="$1"
FASTQ_MANIFEST="$2"
SAMPLES_MANIFEST="$3"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.env
source "$SCRIPT_DIR/config.env"

"$SCRIPT_DIR/download_project.sh" "$PROJECT" "$FASTQ_MANIFEST"

while [[ ! -f "$RNA_ROOT/state/$REFERENCE_COMPLETE_MARKER" ]]; do
  echo "Waiting for reference index: $(date --iso-8601=seconds)"
  sleep 30
done

"$SCRIPT_DIR/quantify_project.sh" "$PROJECT" "$SAMPLES_MANIFEST"
