#!/usr/bin/env bash
set -euo pipefail

# gapless (Nanopore) round starting from the PacBio gapless assembly.
#
# Inputs:
#   gapless_pb/gapless.fa
#   fastq_pass/all_nanopore.fastq.gz
# Outputs:
#   gapless_np/gapless.fa (symlink to pass3/gapless.fa)

source "$(dirname "$0")/_common.sh"

THREADS="${THREADS:-80}"
ASSEMBLY_IN="${ASSEMBLY_IN:-$ROOT_DIR/gapless_pb/gapless.fa}"
ONT_FASTQ="${ONT_FASTQ:-$ROOT_DIR/fastq_pass/all_nanopore.fastq.gz}"
OUTDIR="${OUTDIR:-$ROOT_DIR/gapless_np}"

require gapless.sh

if [[ ! -s "$ASSEMBLY_IN" ]]; then
  echo "[ERROR] Missing assembly: $ASSEMBLY_IN" >&2
  exit 2
fi
if [[ ! -s "$ONT_FASTQ" ]]; then
  echo "[ERROR] Missing ONT FASTQ: $ONT_FASTQ" >&2
  exit 3
fi

rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

log "gapless nanopore -> $OUTDIR/gapless.fa"
gapless.sh -j "$THREADS" -i "$ASSEMBLY_IN" -t nanopore -o "$OUTDIR" "$ONT_FASTQ"

[[ -s "$OUTDIR/gapless.fa" ]] || {
  echo "[ERROR] gapless did not produce $OUTDIR/gapless.fa" >&2
  exit 4
}

log "DONE: $OUTDIR/gapless.fa"
