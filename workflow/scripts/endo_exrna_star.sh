# https://github.com/gersteinlab/exceRpt/blob/e8fe71c42777366e4b2bf8e52854d29b74721b5d/ExampleData/testData_human.fastq/endogenousAlignments_genome_Log.out

genome_dir=$1
fq=$2
threads=$3
output=$4
out_prefix=$(echo $output | sed 's/\.*Aligned.out.bam$//g')

STAR \
    --genomeDir $genome_dir \
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
