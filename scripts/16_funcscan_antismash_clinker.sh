#!/usr/bin/env bash
set -euo pipefail

# Run nf-core/funcscan (antiSMASH + ARG/BGC screening) and generate clinker plots.
#
# Requirements:
#   - nextflow
#   - docker/singularity (depending on nf-core profile)
#   - clinker (pip install clinker)
#
# Notes:
#   This step is optional and can take significant compute/storage.

source "$(dirname "$0")/_common.sh"

THREADS="${THREADS:-16}"
PROFILE="${PROFILE:-docker}"
OUTDIR="${OUTDIR:-$ROOT_DIR/funcscan_results}"
SAMPLESHEET="${SAMPLESHEET:-$ROOT_DIR/samplesheet_for_funcscan.csv}"

require nextflow

log "Running nf-core/funcscan ..."
nextflow run nf-core/funcscan \
  -profile "$PROFILE" \
  --input "$SAMPLESHEET" \
  --outdir "$OUTDIR" \
  --run_bgc_screening \
  --run_arg_screening

# Example clinker usage (edit to match your antiSMASH output path)
# ANTISMASH_DIR="$OUTDIR/bgc/antismash/lferri"
# clinker "$ANTISMASH_DIR"/NZ_*.gbk -p lferri_all.html

log "DONE: $OUTDIR"
