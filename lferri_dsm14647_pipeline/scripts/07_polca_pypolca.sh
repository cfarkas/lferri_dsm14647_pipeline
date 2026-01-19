#!/usr/bin/env bash
set -euo pipefail

# Polish the assembly using PyPolca (POLCA backend) with Illumina paired-end reads.
#
# Inputs:
#   canu_lr/lferri.contigs.fasta (default; override ASSEMBLY)
#   fastq/*R1*.fastq.gz and fastq/*R2*.fastq.gz (polishing reads)
# Outputs:
#   polca_out/pypolca_corrected.fasta

source "$(dirname "$0")/_common.sh"

THREADS="${THREADS:-80}"
ASSEMBLY="${ASSEMBLY:-$ROOT_DIR/canu_lr/lferri.contigs.fasta}"
OUTDIR="${OUTDIR:-$ROOT_DIR/polca_out}"
FASTQ_DIR="${FASTQ_DIR:-$ROOT_DIR/fastq}"

require pypolca

if [[ ! -s "$ASSEMBLY" ]]; then
  echo "[ERROR] Missing assembly: $ASSEMBLY" >&2
  exit 1
fi

mkdir -p "$FASTQ_DIR"

R1_ALL="$FASTQ_DIR/allR1.fq.gz"
R2_ALL="$FASTQ_DIR/allR2.fq.gz"

# Concatenate all polishing read pairs (genomic Illumina)
# Adjust the glob if your file naming differs.
cat "$FASTQ_DIR"/*R1*fastq.gz > "$R1_ALL"
cat "$FASTQ_DIR"/*R2*fastq.gz > "$R2_ALL"

rm -rf "$OUTDIR"

pypolca run \
  -a "$ASSEMBLY" \
  -1 "$R1_ALL" \
  -2 "$R2_ALL" \
  -t "$THREADS" \
  -o "$OUTDIR" \
  --careful --memory_limit 4G --force

FINAL_ASM="$OUTDIR/pypolca_corrected.fasta"
if [[ ! -s "$FINAL_ASM" ]]; then
  echo "[ERROR] POLCA did not produce $FINAL_ASM" >&2
  exit 2
fi

echo "[INFO] Polished assembly: $FINAL_ASM"
ls -lh "$FINAL_ASM"
