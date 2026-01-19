#!/usr/bin/env bash
set -euo pipefail

# gapless (PacBio CLR) round starting from original reference.
#
# Inputs:
#   ref/lferri.fa
#   DSM14647_PacBio_native/fastq/all_pacbio_subreads.fastq.gz
# Outputs:
#   gapless_pb/gapless.fa (symlink to pass3/gapless.fa)

source "$(dirname "$0")/_common.sh"

THREADS="${THREADS:-80}"
ASSEMBLY_IN="${ASSEMBLY_IN:-$ROOT_DIR/ref/lferri.fa}"
PB_FASTQ="${PB_FASTQ:-$ROOT_DIR/DSM14647_PacBio_native/fastq/all_pacbio_subreads.fastq.gz}"
OUTDIR="${OUTDIR:-$ROOT_DIR/gapless_pb}"

require gapless.sh

mkdir -p "$OUTDIR"
rm -rf "$OUTDIR"/pass* "$OUTDIR"/logs "$OUTDIR"/timing 2>/dev/null || true

log "gapless pb_clr -> $OUTDIR"
gapless.sh -j "$THREADS" -i "$ASSEMBLY_IN" -t pb_clr -o "$OUTDIR" "$PB_FASTQ"

log "DONE"
