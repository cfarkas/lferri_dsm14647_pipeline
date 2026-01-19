#!/usr/bin/env bash
set -euo pipefail

# Trim Illumina RNA-seq reads with fastp + run FastQC/MultiQC.
#
# Inputs:
#   samples.tsv (default: samples.tsv at repo root; or config/samples.tsv)
# Outputs:
#   fastq/<sample>_R1.trim.fq.gz
#   fastq/<sample>_R2.trim.fq.gz
#   qc/<sample>_fastp.{html,json}
#   qc/*fastqc.{html,zip}
#   qc/multiqc_report.html

source "$(dirname "$0")/_common.sh"

THREADS="${THREADS:-16}"
SAMPLES_TSV="${SAMPLES_TSV:-$ROOT_DIR/samples.tsv}"
if [[ ! -s "$SAMPLES_TSV" ]]; then
  if [[ -s "$ROOT_DIR/config/samples.tsv" ]]; then
    SAMPLES_TSV="$ROOT_DIR/config/samples.tsv"
  else
    echo "[ERROR] samples.tsv not found. Put it at repo root or config/samples.tsv" >&2
    exit 1
  fi
fi

require fastp
require fastqc
require multiqc

mkdir -p "$ROOT_DIR/fastq" "$ROOT_DIR/qc"

# Skip header
{
  read -r header
  while IFS=$'\t' read -r sample_id condition r1 r2; do
    [[ -z "${sample_id:-}" ]] && continue

    out1="$ROOT_DIR/fastq/${sample_id}_R1.trim.fq.gz"
    out2="$ROOT_DIR/fastq/${sample_id}_R2.trim.fq.gz"

    log "fastp: $sample_id ($condition)"
    fastp \
      -i "$r1" -I "$r2" \
      -o "$out1" -O "$out2" \
      --thread "$THREADS" \
      --html "$ROOT_DIR/qc/${sample_id}_fastp.html" \
      --json "$ROOT_DIR/qc/${sample_id}_fastp.json"

    log "FastQC: $sample_id"
    fastqc -t "$THREADS" -o "$ROOT_DIR/qc" "$out1" "$out2"

  done
} < "$SAMPLES_TSV"

log "MultiQC"
multiqc -o "$ROOT_DIR/qc" "$ROOT_DIR/qc"

log "DONE"
