#!/usr/bin/Rscript

library(biomaRt)
library(optparse)
library(stringr)

option_list = list(
  make_option(c("-g", "--genes"), type="character", default=NULL, 
              help="file with genes of interest", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

opt$genes -> genes
opt$out -> outfile
#genes <- "3532_genes.txt"
genes_in <- read.table(genes, header = F)$V1
sample_name <- str_split_fixed(genes, pattern = "_", n = 2)[,1]

bed_from_vector <- function(input){
  # take in vector of hgnc symbols, return dataframe of chromosome coordinates +- 5kb in bed format with header
    mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="https://www.ensembl.org")
  info <- getBM(attributes=c("hgnc_symbol",
                             "chromosome_name",
                             "start_position",
                             "end_position"
  ),
  filters = "hgnc_symbol",
  values = input,
  mart = mart,
  useCache=FALSE)
outfile <- data.frame("#chr" = info$chromosome_name, "start" = info$start_position-5000, "end" = info$end_position + 5000, check.names = FALSE)
return(outfile)
}
bed_3532 <- bed_from_vector(genes_in  )
outfile_name <- paste(sample_name, ".bed", sep = "")
write.table(bed_3532, file= outfile_name, quote = F, row.names = F)
