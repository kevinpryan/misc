#!/usr/bin/Rscript

######################################################################################
### Apply ssGSEA to GTEx TPM data ###
######################################################################################
# will read in GTEx TPM data along with CTA gene list, summarise duplicated hgnc_symbols in the GTEx list,
# remove sex chromosome genes from the analysis, remove genes with zero expression in >95% of the samples,
# and apply ssGSEA to each sample

# out gives a predefined path/name to the output file. (Default is  "CTA_ssGSEA_out.txt")

# out should contain a header row with the sample names, the first col of the second row is "CTA" 
# and the rest of the cols are ssGSEA outputs for each sample

# example call:
# Rscript 2023-10-10-cta_analysis_gtex_all_samples_ssgsea.R -t gtex_tpm.gct -g CTA_list.txt -o 2023-10-10-CTA_ssGSEA_out.txt


##################################################
# read input arguments

library(optparse)
option_list = list(
  make_option(c("-t", "--tpm"), type="character", default=NULL, 
              help="TPM file", metavar="character"),
  make_option(c("-g", "--genes"), type="character", default=NULL, 
              help="CTA gene set file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="CTA_ssGSEA_out.txt",
              help="specify outfile name", metavar="character")
  
); 
opt <- parse_args(OptionParser(option_list=option_list))
if (length(opt) < 2){
  print_help(opt_parser)
  stop("Must supply at least TPM file and CTA file", call.=FALSE)
}

library(vroom)
library(dplyr)
library(biomaRt)
library(GSVA)

# define filter function
rem <- function(x){
  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))
  r <- as.numeric(apply(x,1,function(i) sum(i == 0) ))
  remove <- which(r > dim(x)[2]*0.95)
  return(remove)
}
# read in files
tpm <- opt$tpm
outfile <- opt$out
genes <- opt$genes
gtex_tpm <- vroom(tpm)#, skip = 2)
# get median for genes that are duplicated
n_occur <- data.frame(table(gtex_tpm$Description))
gtex_tpm <- gtex_tpm %>% dplyr::select(-Name)
dups <- gtex_tpm[gtex_tpm$Description %in% n_occur$Var1[n_occur$Freq > 1],]
not_dups <- gtex_tpm[!(gtex_tpm$Description %in% n_occur$Var1[n_occur$Freq > 1]),]
rm(gtex_tpm)
dups_summarise <- dups %>% group_by(Description) %>% dplyr::summarise((across(where(is.numeric), median)))
dups_removed_tpm <- rbind.data.frame(not_dups, dups_summarise)
rm(dups_summarise)
rm(not_dups)
gc()
# remove X and Y chrom genes from TPM
tpm_genes <- dups_removed_tpm$Description
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl") #, host="uswest.ensembl.org")
info_tpm <- getBM(attributes=c("hgnc_symbol",
                               "ensembl_gene_id",
                               "ensembl_transcript_id",
                               "chromosome_name"),
                  filters = c("hgnc_symbol"),
                  values = tpm_genes,
                  mart = mart,
                  useCache=FALSE,
                  uniqueRows = FALSE) 
info_tpm_remove_sex_chroms_distinct <- info_tpm %>% dplyr::distinct(hgnc_symbol, .keep_all = T) %>%  dplyr::filter(!chromosome_name %in% c("X", "Y"))
dups_removed_tpm_remove_sex_chrom <- dups_removed_tpm %>% dplyr::filter(Description %in% info_tpm_remove_sex_chroms_distinct$hgnc_symbol)
dups_removed_tpm_remove_sex_chrom <- as.data.frame(dups_removed_tpm_remove_sex_chrom)
rownames(dups_removed_tpm_remove_sex_chrom) <- dups_removed_tpm_remove_sex_chrom$Description
dups_removed_tpm_remove_sex_chrom <- dups_removed_tpm_remove_sex_chrom %>% dplyr::select(-Description)
remove <- rem(dups_removed_tpm_remove_sex_chrom)
dups_removed_tpm_remove_sex_chrom <- as.matrix(dups_removed_tpm_remove_sex_chrom[-remove,])

# read in cta genes, get chromosome number and remove X and Y chromosome genes
cta_genes <- unique(read.table(genes)$V1)
info <- getBM(attributes=c("hgnc_symbol",
                           "ensembl_gene_id",
                           "ensembl_transcript_id",
                           "chromosome_name"),
              filters = c("hgnc_symbol"),
              values = cta_genes,
              mart = mart,
              useCache=FALSE,
              uniqueRows = FALSE)
info_remove_sex_chroms <- info %>% dplyr::distinct(hgnc_symbol, .keep_all = T) %>%  dplyr::filter(!chromosome_name %in% c("X", "Y"))
info_with_sex_chroms_distinct <- info %>% dplyr::distinct(hgnc_symbol, .keep_all = T)
cta_genes_list <- list(CTA = info_with_sex_chroms_distinct$hgnc_symbol)

# carry out ssGSEA
cta_es_ssgsea <- gsva(expr = dups_removed_tpm_remove_sex_chrom,
                      gset.idx.list = cta_genes_list,
                      method = "ssgsea"
                      )

write.table(cta_es_ssgsea, file = outfile, quote = F, sep = "\t", row.names = T)
