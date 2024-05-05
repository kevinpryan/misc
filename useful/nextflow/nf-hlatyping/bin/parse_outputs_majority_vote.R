#!/usr/bin/env Rscript

# script to take in the outputs of 4 HLA typing tools and carry out majority voting to get the MHC I calls
# usage: Rscript parse_outputs_majority_vote.R --samplename meta.id --optitype /path/to/optitype/ --polysolver /path/to/polysolver/ --kourami /path/to/kourami/ --hlala /path/to/hlala/ --benchmark /path/to/benchmark/results.txt
# load libraries
library(optparse)
library(dplyr)
library(stringr)
library(tibble)
library(tidyr)
library(vroom)
# source functions
# source("~/Documents/PhD/misc/useful/nextflow/nf-hlatyping/bin/HLA-LA_conversion.R")
# source("~/Documents/PhD/misc/useful/nextflow/nf-hlatyping/bin/Optitype_conversion.R")
# source("~/Documents/PhD/misc/useful/nextflow/nf-hlatyping/bin/Polysolver_conversion.R")
# source("~/Documents/PhD/misc/useful/nextflow/nf-hlatyping/bin/kourami_conversion.R")
# source("~/Documents/PhD/misc/useful/nextflow/nf-hlatyping/bin/majority_voting.R")
# source("~/Documents/PhD/misc/useful/nextflow/nf-hlatyping/bin/df_to_list.R")
# source("~/Documents/PhD/misc/useful/nextflow/nf-hlatyping/bin/are_vectors_identical.R")

source("HLA-LA_conversion.R")
source("Optitype_conversion.R")
source("Polysolver_conversion.R")
source("kourami_conversion.R")
source("majority_voting.R")
source("df_to_list.R")
source("are_vectors_identical.R")
# take in command line arguments

option_list = list(
  make_option(c("-s", "--samplename"), type="character", default=NULL, 
              help="name of sample", metavar="character"),
  make_option(c("-o", "--optitype"), type="character", default=NULL, 
              help="path to optitype output", metavar="character"),
  make_option(c("-p", "--polysolver"), type="character", default=NULL,
              help="path to polysolver output", metavar="character"),
  make_option(c("-l", "--hlala"), type="character", default=NULL,
              help="path to hlala output", metavar="character"),
  make_option(c("-k", "--kourami"), type="character", default=NULL,
              help="path to kourami output", metavar="character"),
  make_option(c("-b", "--benchmark"), type="character", default=NULL,
             help="path to benchmark rankings", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (length(opt) < 6){
  print_help(opt_parser)
  stop("6 arguments must be supplied", call.=FALSE)
}

# assign arguments to variables
opt$samplename -> samplename
opt$optitype -> optitype_in
opt$polysolver -> polysolver_in
opt$hlala -> hlala_in
opt$kourami -> kourami_in
opt$benchmark -> benchmark_in

# Read in files

#samplename <- "3532"
hlala <- toolOutputToR.HLA_LA(hlala_in, mhci_only = T, trim = T)
print(hlala)
optitype <- toolOutputToR.Optitype(optitype_in)
polysolver <- toolOutputToR.Polysolver(polysolver_in, trim = T)
kourami <- toolOutputToR.kourami(kourami_in, mhci_only = T, trim = T)
combined <- rbind(hlala, optitype, polysolver, kourami)
rownames(combined) <- c("hlala", "optitype", "polysolver", "kourami")
combined$tool <- rownames(combined)
combined$sample <- rep(samplename, nrow(combined))

# Read in benchmarking (might not be necessary - just using optitype as best)
benchmark <- read.csv(benchmark_in)
benchmark <- benchmark %>% dplyr::filter(tool %in% c("HLA*LA", "Kourami", "Optitype", "Polysolver") & seq_type == "WES")
benchmark$A <- rev(rank(benchmark$A))
benchmark$B <- rev(rank(benchmark$B))
benchmark$C <- rev(rank(benchmark$C))
benchmark <- benchmark %>% dplyr::select(c(tool, A, B, C))
benchmark$tool <- c("hlala", "kourami", "optitype", "polysolver")

# rename cols
colnames(combined) <- c("A1", "A2", "B1", "B2", "C1", "C2", "tool", "sample")

# Run majority voting for HLA-A
A_list <- df_to_list(combined, cols = c("A1", "A2"))
A_identical <- outer(A_list, A_list, FUN = are_vectors_identical_vectorised)
A_vote <- majority_vote(A_identical, A_list)

# Run majority voting for HLA-B
B_list <- df_to_list(combined, cols = c("B1", "B2"))
B_identical <- outer(B_list, B_list, FUN = are_vectors_identical_vectorised)
B_vote <- majority_vote(B_identical, B_list)

# Run majority voting for HLA-C
C_list <- df_to_list(combined, cols = c("C1", "C2"))
C_identical <- outer(C_list, C_list, FUN = are_vectors_identical_vectorised)
C_vote <- majority_vote(C_identical, C_list)

# prepare outputs: table with all calls across all tools, table with majority voting result
rownames(combined) <- NULL
full_output <- combined %>% relocate(., sample, .before = A1) %>% relocate(., tool, .before = A1)
write.table(full_output, quote = F, row.names = F, sep = "\t", file = paste(samplename, "_all_calls_mhci.tsv", sep = ""))

majority_output <- data.frame(sample = samplename)
majority_output$A1 <- A_vote["A1"]
majority_output$A2 <- A_vote["A2"]
majority_output$B1 <- B_vote["B1"]
majority_output$B2 <- B_vote["B2"]
majority_output$C1 <- C_vote["C1"]
majority_output$C2 <- C_vote["C2"]
write.table(majority_output, quote = F, row.names = F, sep = "\t", file = paste(samplename, "_majority_vote_mhci.tsv", sep = ""))
