Note: STAR-generated files (suffix array, unmapped reads, etc.) are too large to store on github even for a very small working example, so this repository integration testing uses an off-repo directory for data. Can still be built on a 4 core in reasonable time.


# Changelog

-   <span class="timestamp-wrapper"><span class="timestamp">[2022-07-28 Thu] </span></span> v1.1: moved symlink command back to integration testing snakefile and set to simple file path. Made reference genome and gtf into params. Added mem resource limit to flexbar.
-   <span class="timestamp-wrapper"><span class="timestamp">[2022-07-26 Tue] </span></span> Minimum viable build. Makes a deseq object starting from fastqs, reference fasta, ENCODE miRNA gtf, and a library tsv.
-   <span class="timestamp-wrapper"><span class="timestamp">[2022-06-08 Wed] </span></span> Repository started
