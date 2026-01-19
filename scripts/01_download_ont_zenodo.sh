#!/usr/bin/env bash
set -euo pipefail

# Download Oxford Nanopore (ONT) reads from Zenodo.
#
# Inputs:
#   config/ont_zenodo_urls.txt  (one URL per line)
# Outputs:
#   fastq_pass/<downloaded fastq.gz>
#   fastq_pass/all_nanopore.fastq.gz (concatenated)

source "$(dirname "$0")/_common.sh"

URLS="${URLS:-$ROOT_DIR/config/ont_zenodo_urls.txt}"
OUTDIR="$ROOT_DIR/fastq_pass"
THREADS="${THREADS:-8}"

mkdir -p "$OUTDIR"

require wget
require cat

if [[ ! -s "$URLS" ]]; then
  echo "ERROR: URL list not found or empty: $URLS" >&2
  echo "Edit config/ont_zenodo_urls.txt with one Zenodo file URL per line." >&2
  exit 1
fi

log "Downloading ONT reads listed in: $URLS"

# Download each URL
n=0
while IFS= read -r url; do
  [[ -z "$url" ]] && continue
  [[ "$url" =~ ^# ]] && continue
  n=$((n+1))

  # Try to infer a filename (works for Zenodo direct links)
  fname=$(basename "${url%%\?*}")
  if [[ -z "$fname" || "$fname" == "download" ]]; then
    fname="ont_${n}.fastq.gz"
  fi

  log "[$n] wget -> $OUTDIR/$fname"
  wget -c -O "$OUTDIR/$fname" "$url"
done < "$URLS"

if [[ $n -eq 0 ]]; then
  echo "ERROR: No URLs found in $URLS" >&2
  exit 2
fi

log "Concatenating -> $OUTDIR/all_nanopore.fastq.gz"
cat "$OUTDIR"/*.fastq.gz > "$OUTDIR/all_nanopore.fastq.gz"

log "DONE: $OUTDIR/all_nanopore.fastq.gz"
