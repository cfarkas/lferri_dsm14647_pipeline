#!/usr/bin/env bash
set -euo pipefail

# Align Illumina RNA-seq reads (paired-end) with HISAT2 and produce sorted BAMs.
#
# Inputs:
#   ref/lferri_final_hisat2 (HISAT2 index prefix)
#   fastq/<sample>_R1.trim.fq.gz
#   fastq/<sample>_R2.trim.fq.gz
# Outputs:
#   align/<sample>.bam
#   align/<sample>.bam.bai

source "$(dirname "$0")/_common.sh"

THREADS="${THREADS:-32}"
SAMPLES_TSV="${SAMPLES_TSV:-$ROOT_DIR/samples.tsv}"
INDEX_PREFIX="${INDEX_PREFIX:-$ROOT_DIR/ref/lferri_final_hisat2}"
OUTDIR="${OUTDIR:-$ROOT_DIR/align}"
mkdir -p "$OUTDIR"

require hisat2
require samtools

if [[ ! -f "$SAMPLES_TSV" ]]; then
  log "samples.tsv not found: $SAMPLES_TSV"
  log "Tip: copy config/samples.tsv -> samples.tsv and edit paths."
  exit 1
fi

log "Reading samples from: $SAMPLES_TSV"

# Skip header
while IFS=$'\t' read -r sample_id condition r1 r2; do
  [[ "$sample_id" == "sample_id" ]] && continue
  [[ -z "${sample_id}" ]] && continue

  r1_trim="$ROOT_DIR/fastq/${sample_id}_R1.trim.fq.gz"
  r2_trim="$ROOT_DIR/fastq/${sample_id}_R2.trim.fq.gz"

  if [[ ! -s "$r1_trim" || ! -s "$r2_trim" ]]; then
    log "Missing trimmed FASTQ for $sample_id:"
    log "  $r1_trim"
    log "  $r2_trim"
    log "Run: scripts/10a_fastp_trim_and_qc.sh"
    exit 2
  fi

  bam="$OUTDIR/${sample_id}.bam"

  log "Aligning $sample_id ($condition)"
  hisat2 -p "$THREADS" --dta \
    -x "$INDEX_PREFIX" \
    -1 "$r1_trim" \
    -2 "$r2_trim" \
    | samtools sort -@ "$THREADS" -o "$bam" -

  samtools index "$bam"
  log "Wrote: $bam"

done < "$SAMPLES_TSV"

log "DONE"
