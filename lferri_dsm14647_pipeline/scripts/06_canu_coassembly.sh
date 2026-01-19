#!/usr/bin/env bash
set -euo pipefail

# Co-assemble ONT + PacBio CLR reads using Canu.
#
# Inputs:
#   fastq_pass/all_nanopore.fastq.gz
#   DSM14647_PacBio_native/fastq/all_pacbio_subreads.fastq.gz
#   ref/lferri.fa  (for genomeSize estimate)
# Outputs:
#   canu_lr/lferri.contigs.fasta
#
# Notes:
#  - Canu may exit non-zero on some small genomes; check canu_lr/canu.log
#  - If Canu fails, you can try Raven as a fallback (see script).

source "$(dirname "$0")/_common.sh"

THREADS="${THREADS:-80}"
ONT_READS="${ONT_READS:-$ROOT_DIR/fastq_pass/all_nanopore.fastq.gz}"
PB_READS="${PB_READS:-$ROOT_DIR/DSM14647_PacBio_native/fastq/all_pacbio_subreads.fastq.gz}"
REF_ORIG="${REF_ORIG:-$ROOT_DIR/ref/lferri.fa}"
OUTDIR="${OUTDIR:-$ROOT_DIR/canu_lr}"

require canu
require seqkit

if [[ ! -s "$ONT_READS" ]]; then
  echo "[ERROR] Missing ONT reads: $ONT_READS" >&2
  exit 31
fi
if [[ ! -s "$PB_READS" ]]; then
  echo "[ERROR] Missing PacBio reads: $PB_READS" >&2
  exit 32
fi
if [[ ! -s "$REF_ORIG" ]]; then
  echo "[ERROR] Missing original reference: $REF_ORIG" >&2
  exit 33
fi

GENOME_BP=$(seqkit stats -T "$REF_ORIG" | awk 'NR==2{print $5}')
if [[ -z "$GENOME_BP" || "$GENOME_BP" -lt 100000 ]]; then
  echo "[WARN] Could not detect genome size from $REF_ORIG; defaulting to 2.60m"
  GENOME_SIZE="2.60m"
else
  GENOME_SIZE=$(awk -v n="$GENOME_BP" 'BEGIN{printf "%.2fm", n/1e6}')
fi

echo "[INFO] Canu genomeSize=${GENOME_SIZE}"
rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

set +e
canu -p lferri -d "$OUTDIR" \
  genomeSize="${GENOME_SIZE}" \
  -nanopore-raw "$ONT_READS" \
  -pacbio-raw "$PB_READS" \
  useGrid=false \
  maxThreads="$THREADS" \
  minReadLength=1000 \
  stopOnLowCoverage=0 \
  saveReads=true 2>&1 | tee "$OUTDIR/canu.log"
CSTAT=${PIPESTATUS[0]}
set -e

if [[ -s "$OUTDIR/lferri.contigs.fasta" ]]; then
  echo "[INFO] Canu contigs: $OUTDIR/lferri.contigs.fasta"
  exit 0
fi

echo "[WARN] Canu did not produce lferri.contigs.fasta (exit=$CSTAT)."

echo "[INFO] Optional fallback: Raven (requires raven-assembler + minimap2)."
echo "       Example (uncomment):"
echo "       raven --threads $THREADS $ONT_READS $PB_READS > $OUTDIR/raven.fasta"
exit 1
