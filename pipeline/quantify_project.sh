#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 PROJECT PROJECT_SAMPLES_TSV" >&2
  exit 2
fi

PROJECT="$1"
SAMPLES_MANIFEST="$2"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.env
source "$SCRIPT_DIR/config.env"

REFERENCE_ROOT="$RNA_ROOT/reference/$REFERENCE_NAME"
INDEX_DIR="$REFERENCE_ROOT/$SALMON_INDEX_NAME"
RAW_DIR="$RNA_ROOT/raw/$PROJECT"
OUTPUT_ROOT="$RNA_ROOT/quant/$PROJECT"
STATE_DIR="$RNA_ROOT/state"
mkdir -p "$OUTPUT_ROOT" "$STATE_DIR"

[[ -f "$STATE_DIR/$REFERENCE_COMPLETE_MARKER" ]] || {
  echo "Reference index is not complete" >&2
  exit 1
}
[[ -f "$STATE_DIR/$PROJECT.download.complete.txt" ]] || {
  echo "FASTQ download is not complete for $PROJECT" >&2
  exit 1
}
[[ -f "$SAMPLES_MANIFEST" ]] || {
  echo "Missing samples manifest: $SAMPLES_MANIFEST" >&2
  exit 1
}

SUCCESS_LIST="$STATE_DIR/$PROJECT.quant.success.tsv"
FAILURE_LIST="$STATE_DIR/$PROJECT.quant.failure.tsv"
if [[ ! -f "$SUCCESS_LIST" ]]; then
  printf 'completed_at\tproject\tsample_accession\tsample_alias\tlayout\toutput\tmapping_rate\tprocessed_fragments\n' > "$SUCCESS_LIST"
fi
if [[ ! -f "$FAILURE_LIST" ]]; then
  printf 'failed_at\tproject\tsample_accession\tsample_alias\tlayout\toutput\texit_code\n' > "$FAILURE_LIST"
fi

sample_number=0
sample_total="$(awk -F '\t' -v project="$PROJECT" 'NR > 1 && $1 == project {n++} END {print n + 0}' "$SAMPLES_MANIFEST")"

# Tabs are IFS whitespace in Bash, so `read` collapses adjacent tabs and shifts
# columns whenever a TSV field is empty. Select only the fields used below and
# separate them with a non-whitespace control character to preserve empties.
while IFS=$'\x1f' read -r project sample_accession sample_alias library_layout run_accessions; do
  [[ "$project" == "$PROJECT" ]] || continue
  sample_number=$((sample_number + 1))
  output="$OUTPUT_ROOT/$sample_accession"
  success_marker="$output/quant.complete.txt"

  if [[ -f "$success_marker" && -s "$output/quant.sf" ]]; then
    echo "[$sample_number/$sample_total] Already complete: $sample_accession"
    continue
  fi

  if [[ -e "$output" ]]; then
    output="$OUTPUT_ROOT/$sample_accession.attempt-$(date -u +%Y%m%dT%H%M%SZ)"
    echo "Preserving prior incomplete output; new attempt: $output"
  fi
  mkdir -p "$output"

  IFS=';' read -r -a runs <<< "$run_accessions"
  mates1=()
  mates2=()
  reads=()
  for run in "${runs[@]}"; do
    if [[ "$library_layout" == "PAIRED" ]]; then
      mate1="$RAW_DIR/${run}_1.fastq.gz"
      mate2="$RAW_DIR/${run}_2.fastq.gz"
      [[ -f "$mate1" && -f "$mate2" ]] || {
        echo "Missing paired FASTQs for $run" >&2
        exit 1
      }
      mates1+=("$mate1")
      mates2+=("$mate2")
    elif [[ "$library_layout" == "SINGLE" ]]; then
      read="$RAW_DIR/${run}.fastq.gz"
      [[ -f "$read" ]] || {
        echo "Missing single-end FASTQ for $run" >&2
        exit 1
      }
      reads+=("$read")
    else
      echo "Unsupported library layout for $sample_accession: $library_layout" >&2
      exit 1
    fi
  done

  echo "[$sample_number/$sample_total] Quantifying $sample_accession ($sample_alias; $library_layout; ${#runs[@]} runs)"
  set +e
  if [[ "$library_layout" == "PAIRED" ]]; then
    salmon quant \
      --index "$INDEX_DIR" \
      --geneMap "$REFERENCE_ROOT/build/tx2gene.tsv" \
      --libType A \
      --mates1 "${mates1[@]}" \
      --mates2 "${mates2[@]}" \
      --threads "$QUANT_THREADS" \
      --validateMappings \
      --seqBias \
      --gcBias \
      --output "$output"
  else
    salmon quant \
      --index "$INDEX_DIR" \
      --geneMap "$REFERENCE_ROOT/build/tx2gene.tsv" \
      --libType A \
      --unmatedReads "${reads[@]}" \
      --fldMean "$SINGLE_END_FLD_MEAN" \
      --fldSD "$SINGLE_END_FLD_SD" \
      --threads "$QUANT_THREADS" \
      --validateMappings \
      --seqBias \
      --output "$output"
  fi
  exit_code=$?
  set -e

  if (( exit_code != 0 )); then
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$(date --iso-8601=seconds)" "$PROJECT" "$sample_accession" "$sample_alias" \
      "$library_layout" "$output" "$exit_code" >> "$FAILURE_LIST"
    echo "Salmon failed for $sample_accession; output preserved at $output" >&2
    exit "$exit_code"
  fi

  meta="$output/aux_info/meta_info.json"
  [[ -s "$output/quant.sf" && -s "$meta" ]] || {
    echo "Salmon returned success but expected output is missing for $sample_accession" >&2
    exit 1
  }
  mapping_rate="$(jq -r '.percent_mapped // .mapping_rate // "unknown"' "$meta")"
  processed="$(jq -r '.num_processed // "unknown"' "$meta")"
  {
    date --iso-8601=seconds
    printf 'project=%s\n' "$PROJECT"
    printf 'sample_accession=%s\n' "$sample_accession"
    printf 'sample_alias=%s\n' "$sample_alias"
    printf 'layout=%s\n' "$library_layout"
    printf 'run_accessions=%s\n' "$run_accessions"
    printf 'mapping_rate=%s\n' "$mapping_rate"
    printf 'processed_fragments=%s\n' "$processed"
  } > "$output/quant.complete.txt"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$(date --iso-8601=seconds)" "$PROJECT" "$sample_accession" "$sample_alias" \
    "$library_layout" "$output" "$mapping_rate" "$processed" >> "$SUCCESS_LIST"
done < <(
  awk -F '\t' 'NR > 1 {printf "%s%c%s%c%s%c%s%c%s\n", $1, 31, $4, 31, $6, 31, $13, 31, $15}' \
    "$SAMPLES_MANIFEST"
)

completed="$(awk -F '\t' -v project="$PROJECT" 'NR > 1 && $2 == project {seen[$3] = 1} END {print length(seen)}' "$SUCCESS_LIST")"
if [[ "$completed" != "$sample_total" ]]; then
  echo "Only $completed of $sample_total samples have successful quantifications" >&2
  exit 1
fi

{
  date --iso-8601=seconds
  printf 'project=%s\n' "$PROJECT"
  printf 'samples=%s\n' "$sample_total"
  printf 'reference=%s\n' "$REFERENCE_NAME"
  printf 'salmon_version=%s\n' "$(salmon --version)"
  du -sh "$OUTPUT_ROOT"
} > "$STATE_DIR/$PROJECT.quant.complete.txt"

cat "$STATE_DIR/$PROJECT.quant.complete.txt"
