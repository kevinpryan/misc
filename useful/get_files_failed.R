#!/usr/bin/Rscript
library(dplyr)
bams_successful <- list.files("/home/administrator/Documents/nf-core-rnaseq-subpopulation-bams-out/star_salmon/", pattern = "*.bam$")
samplesheet <- read.csv("/mnt/data/rnaseq_caf_subpopulation_samplesheet.csv")
bams_successful_samples <- gsub('.{19}$', '', bams_successful) # strip .markdup.sorted.bam from bam files
samplesheet_rerun <- samplesheet %>% filter(!(sample %in% bams_successful_samples)) # samples to rerun
write.csv(samplesheet_rerun, file = "/home/administrator/Documents/nf_core_rnaseq_subpopulations_samplesheet_rerun_failed_samples_20230515.csv", quote = F, row.names = F)

