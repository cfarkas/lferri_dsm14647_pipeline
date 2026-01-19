#!/usr/bin/env bash
set -euo pipefail

# Liftoff: lift GenBank (e.g., PacBio annotation) to the polished assembly
# and (optionally) map old gene IDs to new Bakta gene IDs + gene_name.
#
# Default inputs:
#   - OLD_GBK: ./LFTS.gbk
#   - NEW_FA : ./polca_out/pypolca_corrected.fasta
#   - NEW_GTF: ./L_ferri_bakta/L_ferriphilum.fixed.gtf
#
# Outputs (default OUTDIR=./L_ferri_bakta/output_liftoff_gbk):
#   - LFTS.ref.fa / LFTS.ref.gff3
#   - liftoff_LFTS_to_pypolca_corrected.gff3
#   - pypolca_corrected_to_LFTS.tsv (mapping TSV used by de_deseq2.R)

source "$(dirname "$0")/_common.sh"

THREADS="${THREADS:-32}"

# ============================================================
# INPUTS (override via env vars)
# ============================================================
OLD_GBK="${OLD_GBK:-$ROOT_DIR/LFTS.gbk}"
NEW_FA="${NEW_FA:-$ROOT_DIR/polca_out/pypolca_corrected.fasta}"
NEW_GTF="${NEW_GTF:-$ROOT_DIR/L_ferri_bakta/L_ferriphilum.fixed.gtf}"

OUTDIR="${OUTDIR:-$ROOT_DIR/L_ferri_bakta/output_liftoff_gbk}"
mkdir -p "$OUTDIR"

# Labels for output names (override by exporting REF_LABEL/TGT_LABEL)
REF_LABEL="${REF_LABEL:-$(basename "$OLD_GBK" | sed 's/\.[^.]*$//')}"
TGT_LABEL="${TGT_LABEL:-$(basename "$NEW_FA" | sed 's/\.[^.]*$//')}"

# ============================================================
# OUTPUTS
# ============================================================
OLD_FA="${OUTDIR}/${REF_LABEL}.ref.fa"
OLD_GFF3="${OUTDIR}/${REF_LABEL}.ref.gff3"
OLD_GENE_MAP_TSV="${OUTDIR}/${REF_LABEL}.gene_id_gene_name.tsv"
OLD_GENE_MAP_JSON="${OUTDIR}/${REF_LABEL}.gene_id_gene_name.json"

LIFTOFF_OUT="${OUTDIR}/liftoff_${REF_LABEL}_to_${TGT_LABEL}.gff3"
LIFTOFF_LOG="${OUTDIR}/liftoff_${REF_LABEL}_to_${TGT_LABEL}.log"

NEW_GENE_MAP_TSV="${OUTDIR}/${TGT_LABEL}.bakta_gene_id_gene_name.tsv"
OLD_ON_NEW_BED="${OUTDIR}/${REF_LABEL}_lifted_CDS.bed"
NEW_BED="${OUTDIR}/${TGT_LABEL}_bakta_CDS.bed"
OVL="${OUTDIR}/${REF_LABEL}_vs_${TGT_LABEL}.overlaps.tsv"

# IMPORTANT:
# This script names maps as "<target>_to_<ref>" and "<ref>_to_<target>".
# For your downstream DE script, you want:
#   pypolca_corrected_to_LFTS.tsv
MAP_TSV="${OUTDIR}/${TGT_LABEL}_to_${REF_LABEL}.tsv"
MAP_TSV_REV="${OUTDIR}/${REF_LABEL}_to_${TGT_LABEL}.tsv"

MAP_JSON_OLD2NAME="${OUTDIR}/${REF_LABEL}_gene_id_gene_name.from_${TGT_LABEL}.json"
MAP_JSON_OLD2NEW="${OUTDIR}/${REF_LABEL}_to_${TGT_LABEL}.json"

# ============================================================
# 0) Dependency checks
# ============================================================
for exe in python liftoff bedtools; do
  require "$exe"
done

python - <<'PY'
import sys
try:
    import Bio  # noqa
except Exception:
    print("ERROR: biopython is not installed in this environment.", file=sys.stderr)
    print("Install with: conda install -c conda-forge biopython", file=sys.stderr)
    raise
PY

# ============================================================
# 1) Convert GenBank -> reference FASTA + GFF3 + OLD gene dictionary
#    - gene_id is taken as locus_tag when present
#    - gene_name prefers /gene (symbol), else /product, else NA
# ============================================================
log "[1/5] Converting GenBank -> FASTA + GFF3 + gene map"
OLD_GBK="$OLD_GBK" OLD_FA="$OLD_FA" OLD_GFF3="$OLD_GFF3" OLD_TSV="$OLD_GENE_MAP_TSV" OLD_JSON="$OLD_GENE_MAP_JSON" \
python - <<'PY'
import os, json
from urllib.parse import quote
from Bio import SeqIO

gbk = os.environ["OLD_GBK"]
fa_out = os.environ["OLD_FA"]
gff_out = os.environ["OLD_GFF3"]
tsv_out = os.environ["OLD_TSV"]
json_out = os.environ["OLD_JSON"]

def q(s: str) -> str:
    s = "" if s is None else str(s)
    return quote(s, safe=":@|,+-._")

def first(quals, key):
    v = quals.get(key)
    if not v:
        return None
    return v[0] if isinstance(v, list) else v

def norm_id(x: str) -> str:
    if not x:
        return ""
    return x.strip().replace(" ", "_")

records = list(SeqIO.parse(gbk, "genbank"))
if not records:
    raise SystemExit(f"ERROR: no GenBank records found in {gbk}")

# FASTA
with open(fa_out, "w") as fa:
    for rec in records:
        fa.write(f">{rec.id}\n")
        seq = str(rec.seq)
        for i in range(0, len(seq), 60):
            fa.write(seq[i:i+60] + "\n")

# collect meta per locus_tag
gene_meta = {}  # gid -> {gene, product, type}

for rec in records:
    for feat in rec.features:
        if feat.type in ("source",):
            continue
        quals = feat.qualifiers
        gid = first(quals, "locus_tag") or first(quals, "gene")
        gid = norm_id(gid) if gid else None
        if not gid:
            continue
        gsym = first(quals, "gene")
        prod = first(quals, "product")
        if gid not in gene_meta:
            gene_meta[gid] = {"gene": None, "product": None, "type": None}
        if gsym and not gene_meta[gid]["gene"]:
            gene_meta[gid]["gene"] = str(gsym).strip()
        if prod and not gene_meta[gid]["product"]:
            gene_meta[gid]["product"] = str(prod).strip()
        if feat.type and not gene_meta[gid]["type"]:
            gene_meta[gid]["type"] = feat.type

# location helpers

def part_coords(loc):
    start = int(loc.start) + 1
    end = int(loc.end)
    strand = loc.strand
    strand_ch = "+" if strand == 1 else "-" if strand == -1 else "."
    return start, end, strand_ch

def span_coords(location):
    parts = list(getattr(location, "parts", [location]))
    starts = [int(p.start) for p in parts]
    ends = [int(p.end) for p in parts]
    strand = getattr(location, "strand", None)
    strand_ch = "+" if strand == 1 else "-" if strand == -1 else "."
    return min(starts) + 1, max(ends), strand_ch, parts

def phase_from_codon_start(quals):
    cs = first(quals, "codon_start")
    try:
        cs = int(cs)
    except Exception:
        cs = 1
    return str((cs - 1) % 3)

# write GFF3
with open(gff_out, "w") as gff:
    gff.write("##gff-version 3\n")
    auto_id = 0

    for rec in records:
        seqid = rec.id
        for feat in rec.features:
            if feat.type in ("source",):
                continue
            ftype = feat.type
            if ftype not in ("gene", "CDS", "tRNA", "rRNA", "ncRNA", "tmRNA"):
                continue

            quals = feat.qualifiers
            gid = first(quals, "locus_tag") or first(quals, "gene")
            gid = norm_id(gid) if gid else None
            if not gid:
                auto_id += 1
                gid = f"GENE_{auto_id:05d}"

            meta = gene_meta.get(gid, {})
            gsym = meta.get("gene") or first(quals, "gene")
            prod = meta.get("product") or first(quals, "product")

            g_start, g_end, g_strand, parts = span_coords(feat.location)

            gene_id = f"gene:{gid}"
            attrs = [f"ID={q(gene_id)}", f"locus_tag={q(gid)}"]
            if gsym:
                attrs.append(f"Name={q(gsym)}")
                attrs.append(f"gene={q(gsym)}")
            if prod:
                attrs.append(f"product={q(prod)}")
            gff.write("\t".join([seqid, "genbank", "gene", str(g_start), str(g_end), ".", g_strand, ".", ";".join(attrs)]) + "\n")

            if ftype == "CDS":
                mrna_id = f"mRNA:{gid}"
                mrna_attrs = [f"ID={q(mrna_id)}", f"Parent={q(gene_id)}", f"locus_tag={q(gid)}"]
                if gsym:
                    mrna_attrs.append(f"Name={q(gsym)}")
                if prod:
                    mrna_attrs.append(f"product={q(prod)}")
                gff.write("\t".join([seqid, "genbank", "mRNA", str(g_start), str(g_end), ".", g_strand, ".", ";".join(mrna_attrs)]) + "\n")

                ph0 = phase_from_codon_start(quals)
                transl_table = first(quals, "transl_table")

                for i, p in enumerate(parts, start=1):
                    p_start, p_end, p_strand = part_coords(p)
                    phase = ph0 if i == 1 else "0"
                    cds_id = f"cds:{gid}:{i}"
                    cds_attrs = [f"ID={q(cds_id)}", f"Parent={q(mrna_id)}", f"locus_tag={q(gid)}"]
                    if gsym:
                        cds_attrs.append(f"gene={q(gsym)}")
                    if prod:
                        cds_attrs.append(f"product={q(prod)}")
                    if transl_table:
                        cds_attrs.append(f"transl_table={q(transl_table)}")
                    gff.write("\t".join([seqid, "genbank", "CDS", str(p_start), str(p_end), ".", p_strand, phase, ";".join(cds_attrs)]) + "\n")

            elif ftype in ("tRNA", "rRNA", "ncRNA", "tmRNA"):
                rna_id = f"{ftype}:{gid}"
                rna_attrs = [f"ID={q(rna_id)}", f"Parent={q(gene_id)}", f"locus_tag={q(gid)}"]
                if gsym:
                    rna_attrs.append(f"Name={q(gsym)}")
                    rna_attrs.append(f"gene={q(gsym)}")
                if prod:
                    rna_attrs.append(f"product={q(prod)}")
                gff.write("\t".join([seqid, "genbank", ftype, str(g_start), str(g_end), ".", g_strand, ".", ";".join(rna_attrs)]) + "\n")

# TSV+JSON dictionary
old_map = {}
with open(tsv_out, "w") as oh:
    oh.write("\t".join(["gene_id", "gene_name", "product"]) + "\n")
    for gid in sorted(gene_meta.keys()):
        gsym = gene_meta[gid].get("gene")
        prod = gene_meta[gid].get("product")
        gname = (gsym.strip() if gsym else None) or (prod.strip() if prod else None) or "NA"
        prod_out = prod.strip() if prod else "NA"
        old_map[gid] = None if gname == "NA" else gname
        oh.write(f"{gid}\t{gname}\t{prod_out}\n")

with open(json_out, "w") as jh:
    json.dump(old_map, jh, indent=2)

print("Wrote:", fa_out)
print("Wrote:", gff_out)
print("Wrote:", tsv_out)
print("Wrote:", json_out)
PY

[[ -s "$OLD_FA" ]] || { echo "ERROR: reference FASTA not created: $OLD_FA" >&2; exit 1; }
[[ -s "$OLD_GFF3" ]] || { echo "ERROR: reference GFF3 not created: $OLD_GFF3" >&2; exit 1; }

# ============================================================
# 2) Liftoff
# ============================================================
log "[2/5] Running Liftoff"
: > "$LIFTOFF_LOG"

run_liftoff() {
  local annot="$1"
  echo "[Liftoff] annotation: $annot" | tee -a "$LIFTOFF_LOG"
  rm -f "$LIFTOFF_OUT" || true

  # Order A
  if liftoff -g "$annot" -o "$LIFTOFF_OUT" -p "$THREADS" "$NEW_FA" "$OLD_FA" >>"$LIFTOFF_LOG" 2>&1; then
    echo "[Liftoff] SUCCESS (order A)" | tee -a "$LIFTOFF_LOG"
    return 0
  fi

  # Order B
  rm -f "$LIFTOFF_OUT" || true
  if liftoff "$NEW_FA" "$OLD_FA" -g "$annot" -o "$LIFTOFF_OUT" -p "$THREADS" >>"$LIFTOFF_LOG" 2>&1; then
    echo "[Liftoff] SUCCESS (order B)" | tee -a "$LIFTOFF_LOG"
    return 0
  fi

  return 1
}

run_liftoff "$OLD_GFF3" || { echo "ERROR: Liftoff failed. See: $LIFTOFF_LOG" >&2; exit 1; }
[[ -s "$LIFTOFF_OUT" ]] || { echo "ERROR: Liftoff output empty: $LIFTOFF_OUT" >&2; exit 1; }

log "Wrote: $LIFTOFF_OUT"
log "Log:   $LIFTOFF_LOG"

# ============================================================
# If NEW_GTF missing: stop here
# ============================================================
if [[ -z "${NEW_GTF:-}" || ! -f "$NEW_GTF" ]]; then
  log "NEW_GTF not provided/found; skipping overlap mapping"
  exit 0
fi

# ============================================================
# 3) Build NEW dictionary: Bakta gene_id -> gene_name
# ============================================================
log "[3/5] Building Bakta gene_id -> gene_name map"

gawk 'BEGIN{FS=OFS="\t"; print "gene_id","gene_name"}
     /^#/ {next}
     ($3=="CDS" || $3=="exon") {
       gid=""; gname="NA";
       if (match($9, /gene_id "([^"]+)"/, m)) gid=m[1]; else next;
       if (match($9, /gene_name "([^"]+)"/, n)) gname=n[1];
       if (!(gid in seen)) {seen[gid]=gname; order[++k]=gid}
       else {if (seen[gid]=="NA" && gname!="NA") seen[gid]=gname}
     }
     END{for(i=1;i<=k;i++){g=order[i]; print g, seen[g]}}' "$NEW_GTF" > "$NEW_GENE_MAP_TSV"

log "Wrote: $NEW_GENE_MAP_TSV"

# ============================================================
# 4) Extract BEDs and overlaps
# ============================================================
log "[4/5] Extracting BEDs + overlaps"

GFF3_IN="$LIFTOFF_OUT" python - <<'PY' > "$OLD_ON_NEW_BED"
import os
from urllib.parse import unquote

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
            label = a.get("locus_tag") or a.get("gene_id") or a.get("Name") or fid
            label = strip_id(label)
            if fid:
                gene_by_id[fid] = label
            continue

        if ftype in {"mRNA", "transcript", "rRNA", "tRNA", "ncRNA"}:
            label = gene_by_id.get(parent, parent) if parent else (a.get("locus_tag") or a.get("gene_id") or a.get("Name") or fid)
            label = strip_id(label)
            if fid:
                tx_by_id[fid] = label
            continue

        if ftype == "CDS":
            gene_label = ""
            if parent:
                gene_label = tx_by_id.get(parent) or gene_by_id.get(parent) or parent
            if not gene_label:
                gene_label = a.get("locus_tag") or a.get("gene_id") or a.get("Name") or parent or fid
            gene_label = strip_id(gene_label)
            if not gene_label:
                continue
            s0 = int(start) - 1
            e1 = int(end)
            print(seqid, s0, e1, gene_label, 0, strand, sep="\t")
PY

NEW_GTF_IN="$NEW_GTF" python - <<'PY' > "$NEW_BED"
import os, re

gtf = os.environ["NEW_GTF_IN"]
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

# Overlap (strand first; fallback unstranded)
bedtools intersect -s -wo -a "$OLD_ON_NEW_BED" -b "$NEW_BED" > "$OVL" || true
if [[ ! -s "$OVL" ]]; then
  log "No overlaps with -s; retrying unstranded"
  bedtools intersect -wo -a "$OLD_ON_NEW_BED" -b "$NEW_BED" > "$OVL"
fi

log "Wrote: $OVL"

# ============================================================
# 5) Pick best matches & write maps
# ============================================================
log "[5/5] Building mapping TSV/JSON"

NEW_GENE_MAP="$NEW_GENE_MAP_TSV" OLD_BED="$OLD_ON_NEW_BED" NEW_BED_IN="$NEW_BED" OVL_IN="$OVL" \
OUT_MAP="$MAP_TSV" OUT_REV="$MAP_TSV_REV" OUT_JSON_NAME="$MAP_JSON_OLD2NAME" OUT_JSON_OLD2NEW="$MAP_JSON_OLD2NEW" \
python - <<'PY'
from collections import defaultdict
import os, json

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

# NEW gene_id -> gene_name
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

# NOTE: out_map is NEW->OLD here (because your downstream wants pypolca_corrected_to_LFTS.tsv)
with open(out_map, "w") as out:
    out.write("\t".join(["new_gene_id","old_gene_id","gene_name","overlap_bp","frac_old","frac_new"]) + "\n")
    for new_id in sorted(new_len.keys()):
        old_id, ovbp = best_old_for_new.get(new_id, ("NA", 0))
        olen = old_len.get(old_id, 0) if old_id != "NA" else 0
        nlen = new_len.get(new_id, 0)

        frac_old = (ovbp / olen) if olen else 0.0
        frac_new = (ovbp / nlen) if nlen else 0.0

        gname = gene_name.get(new_id)

        # Apply minimal overlap filter only for mapping confidence
        if ovbp < MIN_OVERLAP_BP or (olen and frac_old < MIN_FRAC_OLD):
            old_id = "NA"
            ovbp = 0
            frac_old = 0.0
            frac_new = 0.0

        out.write("\t".join([
            new_id,
            old_id,
            (gname if gname is not None else "NA"),
            str(ovbp),
            f"{frac_old:.4f}",
            f"{frac_new:.4f}",
        ]) + "\n")

with open(out_rev, "w") as out:
    out.write("\t".join(["old_gene_id","new_gene_id","gene_name","overlap_bp","frac_old","frac_new"]) + "\n")
    for old_id in sorted(old_len.keys()):
        new_id, ovbp = best_new_for_old.get(old_id, ("NA", 0))
        olen = old_len.get(old_id, 0)
        nlen = new_len.get(new_id, 0) if new_id != "NA" else 0

        frac_old = (ovbp / olen) if olen else 0.0
        frac_new = (ovbp / nlen) if nlen else 0.0

        gname = gene_name.get(new_id)

        if ovbp < MIN_OVERLAP_BP or frac_old < MIN_FRAC_OLD:
            new_id = "NA"
            ovbp = 0
            frac_old = 0.0
            frac_new = 0.0
            gname = None

        out.write("\t".join([
            old_id,
            new_id,
            (gname if gname is not None else "NA"),
            str(ovbp),
            f"{frac_old:.4f}",
            f"{frac_new:.4f}",
        ]) + "\n")

# JSON dictionaries for old gene symbol mapping
old2name = {}
old2new = {}
with open(out_rev) as fh:
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
log "Key mapping TSV for DE script: $MAP_TSV"
