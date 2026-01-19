#!/usr/bin/env bash
set -euo pipefail

# Liftoff: lift OLD RefSeq annotation (from the original NCBI reference) onto
# the NEW polished assembly; then map lifted CDS back to NEW Bakta gene IDs.
#
# This is optional for the RNA-seq pipeline, but useful if you want to translate
# old locus tags (e.g., LPTCAG_RS...) to new Bakta IDs.

source "$(dirname "$0")/_common.sh"

THREADS="${THREADS:-32}"

OLD_GTF="${OLD_GTF:-$ROOT_DIR/ref/lferri.gtf}"
OLD_FA="${OLD_FA:-$ROOT_DIR/ref/lferri.fa}"

NEW_FA="${NEW_FA:-$ROOT_DIR/polca_out/pypolca_corrected.fasta}"
NEW_GTF="${NEW_GTF:-$ROOT_DIR/L_ferri_bakta/L_ferriphilum.fixed.gtf}"

OUTDIR="${OUTDIR:-$ROOT_DIR/L_ferri_bakta/output_liftoff_LPTCAG}"
mkdir -p "$OUTDIR"

require liftoff
require bedtools
require python
require gawk

# -----------------------------
# OUTPUT FILES
# -----------------------------
NEW_GENE_MAP_TSV="${OUTDIR}/gene_id_gene_name.tsv"

LIFTOFF_OUT="${OUTDIR}/liftoff_refseq_to_polca.gff3"
LIFTOFF_LOG="${OUTDIR}/liftoff_refseq_to_polca.log"

OLD_ON_NEW_BED="${OUTDIR}/refseq_LPTCAG_lifted_CDS.bed"
NEW_BED="${OUTDIR}/bakta_IOPOBE_CDS.bed"

OVL="${OUTDIR}/LPTCAG_vs_IOPOBE.overlaps.tsv"

MAP_TSV="${OUTDIR}/LPTCAG_to_IOPOBE_gene_name.tsv"
MAP_TSV_REV="${OUTDIR}/IOPOBE_to_LPTCAG.tsv"

MAP_JSON_OLD2NAME="${OUTDIR}/LPTCAG_gene_id_gene_name.json"
MAP_JSON_OLD2NEW="${OUTDIR}/LPTCAG_to_IOPOBE.json"

# ============================================================
# 1) Build NEW dictionary: IOPOBE gene_id -> gene_name (if present)
# ============================================================
log "Building NEW gene_id -> gene_name map from: $NEW_GTF"

gawk 'BEGIN{FS=OFS="\t"; print "gene_id","gene_name"}
     /^#/ {next}
     ($3=="CDS" || $3=="exon") {
       gid=""; gname="NA";
       if (match($9, /gene_id "([^"]+)"/, m)) gid=m[1]; else next;
       if (match($9, /gene_name "([^"]+)"/, n)) gname=n[1];

       if (!(gid in seen)) {
         seen[gid]=gname;
         order[++k]=gid;
       } else {
         if (seen[gid]=="NA" && gname!="NA") seen[gid]=gname;
       }
     }
     END{for(i=1;i<=k;i++){g=order[i]; print g, seen[g];}}' "$NEW_GTF" > "$NEW_GENE_MAP_TSV"

log "Wrote: $NEW_GENE_MAP_TSV"

# ============================================================
# 2) Liftoff
# ============================================================
: > "$LIFTOFF_LOG"
log "Running Liftoff (RefSeq -> POLCA)"
if ! liftoff -g "$OLD_GTF" -o "$LIFTOFF_OUT" -p "$THREADS" "$NEW_FA" "$OLD_FA" >>"$LIFTOFF_LOG" 2>&1; then
  log "Liftoff failed using GTF directly. Try converting to GFF3 with gffread." 
  require gffread
  OLD_GFF3="${OUTDIR}/lferri.refseq.gff3"
  gffread -E "$OLD_GTF" -o "$OLD_GFF3" >>"$LIFTOFF_LOG" 2>&1
  liftoff -g "$OLD_GFF3" -o "$LIFTOFF_OUT" -p "$THREADS" "$NEW_FA" "$OLD_FA" >>"$LIFTOFF_LOG" 2>&1
fi

[ -s "$LIFTOFF_OUT" ] || { echo "ERROR: Liftoff output is empty: $LIFTOFF_OUT" >&2; exit 1; }
log "Wrote: $LIFTOFF_OUT"

# ============================================================
# 3) Extract BEDs and overlap
# ============================================================
log "Extracting CDS BEDs and computing overlaps"

GFF3_IN="$LIFTOFF_OUT" python - <<'PY' > "$OLD_ON_NEW_BED"
from urllib.parse import unquote
import os

gff3 = os.environ["GFF3_IN"]

def parse_attrs(s: str):
    d = {}
    s = s.strip().strip(";")
    if not s:
        return d
    for part in s.split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = unquote(v)
        elif " " in part:
            k, v = part.split(" ", 1)
            d[k] = v.strip().strip('"')
    return d

def strip_id(x: str) -> str:
    if not x:
        return ""
    x = x.split(",")[0].strip()
    x = unquote(x)
    for p in ("gene:", "transcript:", "mrna:", "mRNA:", "rna:", "cds:", "CDS:"):
        if x.startswith(p):
            x = x[len(p):]
    return x

gene_by_id = {}
tx_by_id = {}

gene_types = {"gene", "pseudogene", "ncRNA_gene", "rRNA_gene", "tRNA_gene"}

with open(gff3) as fh:
    for line in fh:
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        seqid, source, ftype, start, end, score, strand, phase, attr = f
        a = parse_attrs(attr)

        fid = strip_id(a.get("ID", ""))
        parent = strip_id(a.get("Parent", ""))

        if ftype in gene_types:
            label = a.get("gene_id") or a.get("locus_tag") or a.get("Name") or fid
            label = strip_id(label)
            if fid:
                gene_by_id[fid] = label
            continue

        if ftype in {"mRNA", "transcript", "rRNA", "tRNA", "ncRNA"}:
            label = gene_by_id.get(parent, parent) if parent else (a.get("gene_id") or a.get("locus_tag") or a.get("Name") or fid)
            label = strip_id(label)
            if fid:
                tx_by_id[fid] = label
            continue

        if ftype == "CDS":
            gene_label = ""
            if parent:
                gene_label = tx_by_id.get(parent) or gene_by_id.get(parent) or parent
            if not gene_label:
                gene_label = a.get("gene_id") or a.get("locus_tag") or a.get("Name") or parent or fid

            gene_label = strip_id(gene_label)
            if not gene_label:
                continue

            s0 = int(start) - 1
            e1 = int(end)
            print(seqid, s0, e1, gene_label, 0, strand, sep="\t")
PY

GTF_IN="$NEW_GTF" python - <<'PY' > "$NEW_BED"
import os, re

gtf = os.environ["GTF_IN"]
gid_re = re.compile(r'gene_id "([^"]+)"')

with open(gtf) as fh:
    for line in fh:
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        seqid, source, ftype, start, end, score, strand, frame, attr = f
        if ftype != "CDS":
            continue
        m = gid_re.search(attr)
        if not m:
            continue
        gid = m.group(1)
        s0 = int(start) - 1
        e1 = int(end)
        print(seqid, s0, e1, gid, 0, strand, sep="\t")
PY

bedtools intersect -s -wo -a "$OLD_ON_NEW_BED" -b "$NEW_BED" > "$OVL" || true
if [ ! -s "$OVL" ]; then
  log "No overlaps with -s; retrying unstranded"
  bedtools intersect -wo -a "$OLD_ON_NEW_BED" -b "$NEW_BED" > "$OVL"
fi

# ============================================================
# 4) Pick best matches
# ============================================================
log "Selecting best overlaps and writing mapping tables"

NEW_GENE_MAP="$NEW_GENE_MAP_TSV" OLD_BED="$OLD_ON_NEW_BED" NEW_BED_IN="$NEW_BED" OVL_IN="$OVL" \
OUT_MAP="$MAP_TSV" OUT_REV="$MAP_TSV_REV" OUT_JSON_NAME="$MAP_JSON_OLD2NAME" OUT_JSON_OLD2NEW="$MAP_JSON_OLD2NEW" \
python - <<'PY'
from collections import defaultdict
import json, os

new_gene_map = os.environ["NEW_GENE_MAP"]
old_bed = os.environ["OLD_BED"]
new_bed = os.environ["NEW_BED_IN"]
ovl = os.environ["OVL_IN"]

out_map = os.environ["OUT_MAP"]
out_rev = os.environ["OUT_REV"]
out_json_name = os.environ["OUT_JSON_NAME"]
out_json_old2new = os.environ["OUT_JSON_OLD2NEW"]

MIN_OVERLAP_BP = 30
MIN_FRAC_OLD = 0.20

def read_lengths(bed_path):
    lens = defaultdict(int)
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            seq, s, e, gid, score, strand = line.rstrip("\n").split("\t")[:6]
            lens[gid] += int(e) - int(s)
    return dict(lens)

gene_name = {}
with open(new_gene_map) as fh:
    next(fh, None)
    for line in fh:
        if not line.strip():
            continue
        gid, gname = line.rstrip("\n").split("\t", 1)
        gene_name[gid] = None if gname == "NA" else gname

old_len = read_lengths(old_bed)
new_len = read_lengths(new_bed)

pair_ov = defaultdict(int)
with open(ovl) as fh:
    for line in fh:
        if not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        old_id = f[3]
        new_id = f[9]
        overlap_bp = int(f[-1])
        pair_ov[(old_id, new_id)] += overlap_bp

best_new_for_old = {}
for (old_id, new_id), ovbp in pair_ov.items():
    if (old_id not in best_new_for_old) or (ovbp > best_new_for_old[old_id][1]):
        best_new_for_old[old_id] = (new_id, ovbp)

best_old_for_new = {}
for (old_id, new_id), ovbp in pair_ov.items():
    if (new_id not in best_old_for_new) or (ovbp > best_old_for_new[new_id][1]):
        best_old_for_new[new_id] = (old_id, ovbp)

with open(out_map, "w") as out:
    out.write("\t".join(["old_gene_id","new_gene_id","gene_name","overlap_bp","frac_old","frac_new"]) + "\n")
    for old_id in sorted(old_len.keys()):
        new_id, ovbp = best_new_for_old.get(old_id, ("NA", 0))
        olen = old_len.get(old_id, 0)
        nlen = new_len.get(new_id, 0) if new_id != "NA" else 0

        frac_old = (ovbp / olen) if olen else 0.0
        frac_new = (ovbp / nlen) if nlen else 0.0

        if ovbp < MIN_OVERLAP_BP or frac_old < MIN_FRAC_OLD:
            new_id = "NA"
            ovbp = 0
            frac_old = 0.0
            frac_new = 0.0
            gname = None
        else:
            gname = gene_name.get(new_id)

        out.write("\t".join([
            old_id,
            new_id,
            (gname if gname is not None else "NA"),
            str(ovbp),
            f"{frac_old:.4f}",
            f"{frac_new:.4f}",
        ]) + "\n")

with open(out_rev, "w") as out:
    out.write("\t".join(["new_gene_id","old_gene_id","gene_name","overlap_bp","frac_old","frac_new"]) + "\n")
    for new_id in sorted(new_len.keys()):
        old_id, ovbp = best_old_for_new.get(new_id, ("NA", 0))
        olen = old_len.get(old_id, 0) if old_id != "NA" else 0
        nlen = new_len.get(new_id, 0)

        frac_old = (ovbp / olen) if olen else 0.0
        frac_new = (ovbp / nlen) if nlen else 0.0
        gname = gene_name.get(new_id)

        out.write("\t".join([
            new_id,
            old_id,
            (gname if gname is not None else "NA"),
            str(ovbp),
            f"{frac_old:.4f}",
            f"{frac_new:.4f}",
        ]) + "\n")

old2name = {}
old2new = {}
with open(out_map) as fh:
    next(fh)
    for line in fh:
        old_id, new_id, gname, ovbp, frac_old, frac_new = line.rstrip("\n").split("\t")
        old2new[old_id] = (None if new_id == "NA" else new_id)
        old2name[old_id] = (None if gname == "NA" else gname)

with open(out_json_name, "w") as oh:
    json.dump(old2name, oh, indent=2)
with open(out_json_old2new, "w") as oh:
    json.dump(old2new, oh, indent=2)

print("Wrote:", out_map)
print("Wrote:", out_rev)
print("Wrote:", out_json_name)
print("Wrote:", out_json_old2new)
PY

log "DONE"
log "Outputs:"
log "  - $LIFTOFF_OUT"
log "  - $MAP_TSV"
