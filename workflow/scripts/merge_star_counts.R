#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)
count_files_string = args[1]
counts_output_tsv = args[2]

library(plyr)
library(tidyverse)

count_files_list = unlist(strsplit(count_files_string, " "))

count_files = lapply(count_files_list, function(x){read_tsv(x)})
names(count_files) = substr(gsub("^.*lib", "lib", count_files_list), 1, 6)

counts = plyr::join_all(count_files, type="full", by = "ensembl")
row.names(counts) = counts$ensembl
counts = counts[,-1]

write.table(counts, file = counts_output_tsv, row.names = TRUE, sep = '\t', quote = F)
