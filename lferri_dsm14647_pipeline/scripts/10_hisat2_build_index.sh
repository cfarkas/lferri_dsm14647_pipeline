#!/usr/bin/env bash
set -euo pipefail

# Build HISAT2 index on the final polished assembly.
#
# Inputs:
#   polca_out/pypolca_corrected.fasta
# Outputs:
#   ref/lferri_final_hisat2.*.ht2

source "$(dirname "$0")/_common.sh"

THREADS="${THREADS:-32}"
FA="${FA:-$ROOT_DIR/polca_out/pypolca_corrected.fasta}"
INDEX_PREFIX="${INDEX_PREFIX:-$ROOT_DIR/ref/lferri_final_hisat2}"

require hisat2-build

log "Building HISAT2 index: $INDEX_PREFIX"
hisat2-build -p "$THREADS" "$FA" "$INDEX_PREFIX"

log "DONE"
