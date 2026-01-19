SHELL := /usr/bin/env bash

# Default threads (override: make THREADS=80 <target>)
THREADS ?= 32

.PHONY: help
help:
	@echo "Targets:"
	@echo "  fetch_ref          - download original reference to ref/"
	@echo "  ont_zenodo          - download ONT reads from Zenodo URLs list"
	@echo "  pacbio_ena          - download PacBio native reads from ENA"
	@echo "  pacbio_fastq        - convert PacBio native to FASTQ"
	@echo "  canu                - Canu co-assembly (ONT+PacBio)"
	@echo "  polca               - PyPolca/POLCA polishing"
	@echo "  bakta               - Bakta annotation"
	@echo "  gtf                 - create fixed GTF for featureCounts"
	@echo "  trim_qc             - fastp + fastqc + multiqc on RNA-seq reads"
	@echo "  hisat2_index         - build HISAT2 index on polished assembly"
	@echo "  align               - align RNA-seq reads to genome"
	@echo "  counts              - featureCounts + clean matrix"
	@echo "  liftoff_gbk          - Liftoff GenBank annotation to new genome"
	@echo "  deseq2              - DESeq2 + ComplexHeatmap (union heatmap)"
	@echo "  funcscan            - nf-core/funcscan + clinker (optional)"

fetch_ref:
	THREADS=$(THREADS) scripts/00_fetch_reference.sh

ont_zenodo:
	THREADS=$(THREADS) scripts/01_download_ont_zenodo.sh

pacbio_ena:
	THREADS=$(THREADS) scripts/02_download_pacbio_native_ena.sh

pacbio_fastq:
	THREADS=$(THREADS) scripts/03_convert_pacbio_native_to_fastq.sh

canu:
	THREADS=$(THREADS) scripts/06_canu_coassembly.sh

polca:
	THREADS=$(THREADS) scripts/07_polca_pypolca.sh

bakta:
	THREADS=$(THREADS) scripts/08_bakta_annotation.sh

gtf:
	scripts/09_format_gtf.sh

trim_qc:
	THREADS=$(THREADS) scripts/10a_fastp_trim_and_qc.sh

hisat2_index:
	THREADS=$(THREADS) scripts/10_hisat2_build_index.sh

align:
	THREADS=$(THREADS) scripts/11_align_rnaseq.sh

counts:
	THREADS=$(THREADS) scripts/12_featurecounts.sh

liftoff_gbk:
	THREADS=$(THREADS) scripts/13_liftoff_from_gbk.sh

deseq2:
	ROOT_DIR=$$(pwd) scripts/15_run_deseq2.sh

funcscan:
	scripts/16_funcscan_antismash_clinker.sh
