# Data sources

This repository **does not ship raw reads**. Instead, it provides scripts + templates to download them (or you can copy them into the expected folders).

## ONT reads (Zenodo)

1. Put direct download URLs (one per line) in: `config/ont_zenodo_urls.txt`
2. Run:

```bash
bash scripts/01_download_ont_zenodo.sh
```

Outputs:
- downloaded files in `fastq_pass/`
- concatenated file: `fastq_pass/all_nanopore.fastq.gz`

## PacBio native subreads (ENA)

The script `scripts/02_download_pacbio_native_ena.sh` queries ENA for run accessions and downloads the native `.bax.h5` files.

Then convert to FASTQ:

```bash
bash scripts/03_convert_pacbio_native_to_fastq.sh
```

Outputs:
- `DSM14647_PacBio_native/fastq/all_pacbio_subreads.fastq.gz`

## Illumina reads

Place Illumina reads under `fastq/`:

- RNA-seq (paired-end): used by `scripts/10a_fastp_trim_and_qc.sh` → `scripts/11_align_rnaseq.sh`
- Polishing reads (paired-end genomic Illumina): used by `scripts/07_polca_pypolca.sh`

If your naming differs, either:
- edit `samples.tsv` / scripts, or
- symlink your files to the expected names.
