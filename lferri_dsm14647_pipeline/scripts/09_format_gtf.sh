#!/usr/bin/env bash
set -euo pipefail

# Convert Bakta GFF3 -> cleaned/fixed GTF for featureCounts.
#
# Outputs:
#   L_ferri_bakta/L_ferriphilum.fixed.gtf
#
# Notes:
#   - Optionally renames contig_N -> tig0000000N. Disable if your BAM uses different names.

source "$(dirname "$0")/_common.sh"

GFF3="${GFF3:-$ROOT_DIR/L_ferri_bakta/L_ferriphilum.gff3}"
OUTDIR="${OUTDIR:-$ROOT_DIR/L_ferri_bakta}"
PREFIX="${PREFIX:-L_ferriphilum}"
RENAME_CONTIGS="${RENAME_CONTIGS:-0}"

require gffread
require awk
require gawk

clean="$OUTDIR/${PREFIX}.clean.gff3"
gtf_raw="$OUTDIR/${PREFIX}.raw.gtf"
gtf_fixed="$OUTDIR/${PREFIX}.fixed.gtf"

if [[ ! -s "$GFF3" ]]; then
  echo "[ERROR] Missing Bakta GFF3: $GFF3" >&2
  exit 1
fi

# 1) Clean unknown strand ? -> . and stop before ##FASTA
awk 'BEGIN{FS=OFS="\t"}
     /^##FASTA/ {exit}
     /^#/ {print; next}
     { if (NF>=7 && $7=="?") $7="."; print }' \
  "$GFF3" > "$clean"

# 2) Convert GFF3 -> GTF
# -E keeps some consistency checks
# -T makes GTF output
# -o output
#
# NOTE: Bakta's GFF3 already contains gene IDs/locus tags in attributes.
gffread -E "$clean" -T -o "$gtf_raw"

# 3) Fix: ensure gene_id and locus_tag exist (from transcript_id)
# Optionally rename contig_N -> tig%08d
if [[ "$RENAME_CONTIGS" == "1" ]]; then
  gawk 'BEGIN{FS=OFS="\t"}
       /^#/ {print; next}
       {
         if ($1 ~ /^contig_[0-9]+$/) {
           n = substr($1, 8)
           $1 = sprintf("tig%08d", n)
         }
         id=""
         if (match($9, /transcript_id "([^"]+)"/, m)) id=m[1]
         if (id != "") {
           if ($9 !~ /gene_id "/)   $9 = $9 " gene_id \"" id "\";"
           if ($9 !~ /locus_tag "/) $9 = $9 " locus_tag \"" id "\";"
         }
         print
       }' "$gtf_raw" > "$gtf_fixed"
else
  gawk 'BEGIN{FS=OFS="\t"}
       /^#/ {print; next}
       {
         id=""
         if (match($9, /transcript_id "([^"]+)"/, m)) id=m[1]
         if (id != "") {
           if ($9 !~ /gene_id "/)   $9 = $9 " gene_id \"" id "\";"
           if ($9 !~ /locus_tag "/) $9 = $9 " locus_tag \"" id "\";"
         }
         print
       }' "$gtf_raw" > "$gtf_fixed"
fi

echo "Wrote: $gtf_fixed"
ls -lh "$gtf_fixed"
