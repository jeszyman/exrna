#!/usr/bin/env bash

# For unit testing
# singularity shell --bind /mnt ~/sing_containers/biotools.1.0.2.sif
# source ./config/bash_src
# mkdir -p ${data_dir}/analysis/exrna/fastqs
# mkdir -p ${data_dir}/logs
# ln -sf ${data_dir}/inputs/10b.H5KGMDRXY_GGCTAC_L002_R1.fastq.gz ${data_dir}/analysis/exrna/fastqs/lib001_raw_R1.fastq.gz
# ln -sf ${data_dir}/inputs/18b.H5KGMDRXY_GTGAAA_L002_R1.fastq.gz ${data_dir}/analysis/exrna/fastqs/lib002_raw_R2.fastq.gz

# in_fastq=${data_dir}/analysis/exrna/fastqs/lib001_raw_R1.fastq.gz
# out_fastq=${data_dir}/analysis/exrna/fastqs/lib001_proc_R1.fastq.gz
# out_fail=${data_dir}/analysis/exrna/fastqs/lib001_fail_R1.fastq.gz
# log_html=${data_dir}/logs/lib001_fastp.html
# log_json=${data_dir}/logs/lib001_fastp.json
# threads=8

# Command line inputs
in_fastq="${1}"
out_fastq="${2}"
out_fail="${3}"
log_html="${4}"
log_json="${5}"
threads="${6}"


# Fastp
# Quality filters are default, but verbose here
# Length limit added appropriate for these small RNA-seq reads
fastp \
    --adapter_sequence "TGGAATTCTCGGGTGCCAAGG" \
    --dont_eval_duplication \
    --failed_out $out_fail \
    --html $log_html \
    --in1 $in_fastq \
    --json $log_json \
    --max_len1 50 \
    --out1 $out_fastq \
    --overrepresentation_analysis \
    --qualified_quality_phred 15 \
    --thread $threads \
    --trim_poly_g \
    --unqualified_percent_limit 40
