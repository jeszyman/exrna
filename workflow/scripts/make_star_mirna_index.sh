#!/usr/bin/env bash

fasta="$1"
gtf="$2"
outdir="$3"
threads="$4"

STAR \
    --runThreadN $threads \
    --runMode genomeGenerate \
    --genomeDir $outdir \
    --sjdbGTFfile $gtf \
    --sjdbOverhang 1 \
    --genomeSAindexNbases 10 \
    --genomeFastaFiles $fasta
