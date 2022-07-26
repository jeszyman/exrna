#!/usr/bin/env bash

index=$1
fq=$2
threads=$3
output=$4
out_prefix=$(echo $output | sed 's/_Reads.*$/_/g')

STAR \
    --genomeDir $index \
    --outFileNamePrefix $out_prefix \
    --outReadsUnmapped Fastx \
    --quantMode GeneCounts \
    --readFilesIn $fq \
    --runThreadN $threads \
    --readFilesCommand zcat \
    --outSAMtype BAM   Unsorted \
    --outSAMattributes All \
    --outSAMunmapped None \
    --outFilterMultimapNmax 1000000 \
    --outFilterMatchNmin 18 \
    --outFilterMatchNminOverLread 0.9 \
    --outFilterMismatchNmax 1 \
    --outFilterMismatchNoverLmax 0.3 \
    --alignIntronMin 2 \
    --alignIntronMax 1 \
    --alignEndsType Local
