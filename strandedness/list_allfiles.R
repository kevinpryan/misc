#!/usr/bin/Rscript

read1 <- list.files(path = "/data/kryan/rna_seq_bc/caf_subtypes/EGAD00001005744", full.names=T, recursive = T, pattern = "R1.fastq.gz$")
read2 <- list.files(path = "/data/kryan/rna_seq_bc/caf_subtypes/EGAD00001005744", full.names=T, recursive = T, pattern = "R2.fastq.gz$")
all_files <- cbind.data.frame(read1, read2)
write.table(all_files, file = "/data/kryan/rna_seq_bc/caf_subtypes/EGAD00001005744/all_files_combined_EGAD00001005744.txt", quote = F, sep = "\t", row.names = F)
