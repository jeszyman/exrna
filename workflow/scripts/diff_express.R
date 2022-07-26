#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
counts_input = args[1]
coldata_input = args[2]
design = args[3]
deseq_output = args[4]

design = "~ cohort"
library(DESeq2)

counts = read.table(counts_input, header = TRUE)

coldata = read.table(coldata_input, header = TRUE)

design = formula(design)

nofilt_dds = DESeqDataSetFromMatrix(countData = counts,
                             colData = coldata,
                             design = design)

# Filter to rows where all samples have counts
keep = rowSums(counts(nofilt_dds)) >= ncol(counts(nofilt_dds))
dds = nofilt_dds[keep,]

dds = DESeq(dds)
nofilt_dds = DESeq(nofilt_dds)

# rlog
rld = rlog(dds)

# Save
save(nofilt_dds,
     dds,
     rld,
     file = deseq_output)
