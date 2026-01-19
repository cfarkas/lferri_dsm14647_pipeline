#!/usr/bin/env bash
set -euo pipefail

# Convert PacBio native subreads (bax.h5) into gzipped FASTQ.
#
# Expected input layout:
#   DSM14647_PacBio_native/*.bax.h5
# Outputs:
#   DSM14647_PacBio_native/bam/*.subreads.bam
#   DSM14647_PacBio_native/fastq/<movie>.fastq.gz
#   DSM14647_PacBio_native/fastq/all_pacbio_subreads.fastq.gz

source "$(dirname "$0")/_common.sh"

INDIR="${INDIR:-$ROOT_DIR/DSM14647_PacBio_native}"
BAMDIR="$INDIR/bam"
FQDIR="$INDIR/fastq"
mkdir -p "$BAMDIR" "$FQDIR"
cd "$INDIR"

require bax2bam
require bam2fastq

# Find movies (group bax files by common prefix up to '_s1_p0')
mapfile -t prefixes < <(
  ls -1 *.bax.h5 2>/dev/null \
  | sed -E 's/_p0\.[0-9]+\.bax\.h5$//' \
  | sort -u
)

if [[ ${#prefixes[@]} -eq 0 ]]; then
  echo "ERROR: No *.bax.h5 found in $INDIR" >&2
  exit 1
fi

log "Detected ${#prefixes[@]} PacBio movies"

for p in "${prefixes[@]}"; do
  log "bax2bam -> $p"
  # collect bax files for this movie
  mapfile -t baxes < <(ls -1 "${p}"_p0.*.bax.h5 2>/dev/null || true)
  if [[ ${#baxes[@]} -eq 0 ]]; then
    echo "WARN: No bax parts found for prefix $p (skipping)" >&2
    continue
  fi

  outprefix="$BAMDIR/$(basename "$p")"
  bax2bam -o "$outprefix" "${baxes[@]}"

  # Convert to FASTQ
  subreads_bam="${outprefix}.subreads.bam"
  if [[ ! -s "$subreads_bam" ]]; then
    echo "ERROR: expected $subreads_bam" >&2
    exit 2
  fi

  base=$(basename "$p")
  bam2fastq --num-threads "$THREADS" -o "$FQDIR/$base" "$subreads_bam"

  if [[ ! -s "$FQDIR/${base}.fastq.gz" ]]; then
    echo "ERROR: expected $FQDIR/${base}.fastq.gz" >&2
    exit 3
  fi

done

log "Concatenating -> $FQDIR/all_pacbio_subreads.fastq.gz"
cat "$FQDIR"/*.fastq.gz > "$FQDIR/all_pacbio_subreads.fastq.gz"

log "DONE"
