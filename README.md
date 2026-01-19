# Leptospirillum ferriphilum DSM 14647 genome polishing + RNA-seq DE pipeline

This repository contains a **reproducible (scripted) workflow** to:

1. Fetch the **original NCBI reference** for Leptospirillum ferriphilum (GCF_000755505.1).
2. Build a **hybrid long-read assembly** (ONT + PacBio CLR) with **Canu**.
3. **Polish** the assembly with Illumina PE using **PyPolca / POLCA**.
4. Annotate the polished genome with **Bakta**.
5. Align RNA-seq reads with **HISAT2**, quantify with **featureCounts**.
6. Run **DESeq2** differential expression (with **apeglm** shrinkage) and generate:
   - MA plots, volcano plots, PCA, dispersions
   - **ComplexHeatmap** plots (Top 50 per contrast + union heatmap across all conditions)
   - CSV exports including **cluster labels** and per-condition means
7. Lift annotations onto the polished genome using **Liftoff** (optional but recommended for gene symbols):
   - From the original RefSeq GTF
   - From a PacBio GenBank (e.g., `LFTS.gbk`)
8. Optional: run **nf-core/funcscan** + **clinker** for BGC/ARG screening and antiSMASH visualization.

The scripts are designed to match (and help reproduce) the flat, analysis-friendly folder layout you already use (e.g. `align/`, `counts/`, `de/`, `L_ferri_bakta/`, `polca_out/`, `fastq_pass/` etc.).

Large inputs/outputs (FASTQs, BAMs, Bakta DB, etc.) are **not tracked in git** (see `.gitignore`).

---

## Expected folder layout

This repo expects (and will create) the following top-level directories:

- `ref/` – reference FASTA/GTF + indexes
- `fastq/` – Illumina reads (RNA-seq + polishing)
- `fastq_pass/` – ONT reads (from Zenodo) + `all_nanopore.fastq.gz`
- `DSM14647_PacBio_native/` – PacBio native downloads + converted FASTQ
- `canu_lr/` – Canu assembly outputs
- `polca_out/` – PyPolca polishing outputs
- `L_ferri_bakta/` – Bakta annotation + Liftoff outputs
- `qc/`, `align/`, `counts/` – RNA-seq QC, BAMs, counts
- `de/` – DE results + figures
- `funcscan_results/` – optional nf-core/funcscan outputs

---

## Quick start

### 1) Clone

```bash
git clone [<this-repo>](https://github.com/cfarkas/lferri_dsm14647_pipeline.git)
cd lferri_dsm14647_pipeline
```

### 2) Configure sample sheet

Edit `samples.tsv` (or `config/samples.tsv`) to point to your RNA-seq FASTQs.

### 3) Download reads

#### ONT (Zenodo, via wget)

Put the Zenodo **direct download URLs** (one per line) into:

- `config/ont_zenodo_urls.txt`

Then run:

```bash
bash scripts/01_download_ont_zenodo.sh
```

This writes downloads to `fastq_pass/` and creates `fastq_pass/all_nanopore.fastq.gz`.

#### PacBio native (ENA)

```bash
bash scripts/02_download_pacbio_native_ena.sh
bash scripts/03_convert_pacbio_native_to_fastq.sh
```

### 4) Assemble + polish + annotate

```bash
bash scripts/00_fetch_reference.sh
bash scripts/06_canu_coassembly.sh
bash scripts/07_polca_pypolca.sh
bash scripts/08_bakta_annotation.sh
bash scripts/09_format_gtf.sh
```

### 5) RNA-seq alignment + counts

```bash
bash scripts/10a_fastp_trim_and_qc.sh
bash scripts/10_hisat2_build_index.sh
bash scripts/11_align_rnaseq.sh
bash scripts/12_featurecounts.sh
```

### 6) Liftoff mapping for gene symbols (recommended)

If you have a GenBank file with gene symbols/products (example: `LFTS.gbk`) in the repo root:

```bash
bash scripts/13_liftoff_from_gbk.sh
```

This produces:

- `L_ferri_bakta/output_liftoff_gbk/pypolca_corrected_to_LFTS.tsv`

which is used by the DE script to append `gene_symbol` to IOPOBE gene IDs.

### 7) DESeq2 + ComplexHeatmap

```bash
bash scripts/15_run_deseq2.sh
```

---

## Outputs

Key outputs:

- `de/DE_*` (tables + significant lists)
- `de/figures/ComplexHeatmap_DEG_union_all_conditions.pdf`
- `de/DEG_union_cluster_condition_means.csv`
- `de/DEG_union_clusters_vst_per_sample.csv`
- `de/DEG_union_clusters_zscore_per_sample.csv`

---

## Citation

If you use this workflow in a publication, please cite the main methods:

- DESeq2 (Love et al., Genome Biology 2014)
- apeglm (Zhu et al., Bioinformatics 2019)
- ComplexHeatmap (Gu et al., Bioinformatics 2016)

And cite the tools you use (HISAT2, subread/featureCounts, Canu, PyPolca/POLCA, Bakta, Liftoff, etc.).

See `docs/references.md` for links.

---

## License

MIT – see `LICENSE`.
