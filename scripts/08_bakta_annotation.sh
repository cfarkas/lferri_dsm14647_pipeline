#!/usr/bin/env bash
set -euo pipefail

# Run Bakta annotation on the polished assembly.
#
# Inputs:
#   polca_out/pypolca_corrected.fasta
# Outputs:
#   L_ferri_bakta/*
#
# Requirements:
#   - bakta installed (and BAKTA_DB env var set to your DB path)
#   - amrfinderplus optional (for update)

source "$(dirname "$0")/_common.sh"

THREADS="${THREADS:-80}"
ASSEMBLY="${ASSEMBLY:-$ROOT_DIR/polca_out/pypolca_corrected.fasta}"
OUTDIR="${OUTDIR:-$ROOT_DIR/L_ferri_bakta}"
PREFIX="${PREFIX:-L_ferriphilum}"

require bakta

if [[ ! -s "$ASSEMBLY" ]]; then
  echo "[ERROR] Missing assembly: $ASSEMBLY" >&2
  exit 1
fi

if [[ -z "${BAKTA_DB:-}" ]]; then
  echo "[ERROR] Please export BAKTA_DB=/path/to/bakta_db before running." >&2
  exit 2
fi

echo "[INFO] Using BAKTA_DB=$BAKTA_DB"

mkdir -p "$OUTDIR"
rm -rf "$OUTDIR"/*

# Optional: update AMRFinder DB (safe to skip)
if command -v amrfinder_update >/dev/null 2>&1; then
  echo "[INFO] Updating AMRFinder database (optional)"
  amrfinder_update --force_update --database "$BAKTA_DB/amrfinderplus-db" || true
fi

bakta --db "$BAKTA_DB" \
  --threads "$THREADS" \
  --prefix "$PREFIX" \
  --output "$OUTDIR" \
  --force \
  "$ASSEMBLY"

echo "[INFO] Bakta output written to: $OUTDIR"
ls -lh "$OUTDIR" | head
