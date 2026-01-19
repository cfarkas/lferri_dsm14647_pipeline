#!/usr/bin/env bash
set -euo pipefail

# Run DESeq2 + ComplexHeatmap pipeline.
#
# Inputs:
#   counts/gene_counts.matrix.tsv
#   samples.tsv
#   L_ferri_bakta/output_liftoff_gbk/pypolca_corrected_to_LFTS.tsv (optional but recommended)
# Outputs:
#   de/*

source "$(dirname "$0")/_common.sh"

require Rscript

ROOT_DIR="${ROOT_DIR:-$ROOT_DIR}"
export ROOT_DIR

Rscript "$ROOT_DIR/de_deseq2.R"
