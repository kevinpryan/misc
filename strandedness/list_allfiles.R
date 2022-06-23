#!/usr/bin/Rscript
library(stringr)
read1 <- list.files(path = "/data/kryan/rna_seq_bc/caf_subtypes/EGAD00001003808", full.names=T, recursive = T, pattern = "R1.fastq.gz$")
read2 <- list.files(path = "/data/kryan/rna_seq_bc/caf_subtypes/EGAD00001003808", full.names=T, recursive = T, pattern = "R2.fastq.gz$")
sample_files <- str_split_fixed(string = read1, pattern = "/", n = 8)[,8]
sample_name <- str_split_fixed(sample_files, pattern = "\\.", n = 4)[,1]
all_files <- cbind.data.frame(sample_name, read1, read2)
write.table(all_files, file = "/data/kryan/rna_seq_bc/caf_subtypes/EGAD00001003808/all_files_combined_EGAD00001003808.txt", quote = F, sep = "\t", row.names = F)
