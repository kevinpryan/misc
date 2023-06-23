#!/usr/bin/Rscript

## Compare CAFs and TANs with GTEx fibroblasts

# read in counts in-house
# read in counts gtex
# match by genes
# make into one deseq object
# batch correction combat-seq
# PCA

### Load libraries
library(tximeta)
library(DESeq2)
library(vroom)
library(sva)
library(dplyr)
library(stringr)
library(BatchQC)
library(optparse)

option_list = list(
  make_option(c("-m", "--method"), type="character", default=NULL, 
              help="nextflow", metavar="character"),
  make_option(c("-t", "--tx2gene_file"), type="character", default=NULL, 
              help="tx2gene file", metavar="character"),
  make_option(c("-i", "--inhouse_metadata"), type="character", default=NULL,
              help="", metavar="character"),
  make_option(c("-d", "--metadata"), type="character", default=NULL,
              help="metadata without pt", metavar="character"),
  make_option(c("-g", "--gtex_file"), type="character", default="eqtl_out.txt",
              help="gtex file", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
params = parse_args(opt_parser);

## Convert command-line arguments to variables
if (params$method == "nextflow") {
  if (length(params) < 4){
    print_help(opt_parser)
    stop("At least 4 arguments must be supplied", call.=FALSE)
  }
  print("params:")
  print(params)
  params$method -> method
  params$tx2gene_file -> tx2gene_file
  params$inhouse_metadata -> inhouse_metadata
  params$metadata -> metadata
  params$gtex_file -> gtex_file
  metadata_original <- read.table(metadata)
  #inhouse_metadata_in <- read.csv(inhouse_metadata)
  metadata_original <- read.table(metadata)
  metadata_original_inhouse <- metadata_original[metadata_original$Study == "InHouse",]
  metadata_original_inhouse$Sample <- as.character(rownames(metadata_original_inhouse))
} else {
  metadata <- "/home/kevin/Documents/PhD/subtypes/caf-subtype-analysis/nf-subpop/outdir/metadata/metadata_full.txt"
  #metadata_path <- "~/Documents/PhD/subtypes/caf-subtype-analysis/nf-subpop/outdir/metadata_by_patient/metadata_with_patient.txt"
  tx2gene_file <- "~/Documents/PhD/subtypes/caf-subtype-analysis/nf-subpop/outdir/tx2gene/tx2gene.txt"
  inhouse_metadata <- "~/Documents/PhD/CAF_data/InHouse/reformat_samples.csv"
  gtex_file <- "~/Downloads/gene_reads_2017-06-05_v8_cells_cultured_fibroblasts.gct/gene_reads_cells_cultured_fibroblasts.gct"
  sub <- "rstudio"

  metadata_original <- read.table(metadata)
  metadata_original_inhouse <- metadata_original[metadata_original$Study == "InHouse",]
  metadata_original_inhouse$Sample <- as.character(rownames(metadata_original_inhouse))
  metadata_original_inhouse$directory <- gsub("kevin", sub, metadata_original_inhouse$directory)
}

inhouse_metadata_in <- read.csv(inhouse_metadata)
colnames(inhouse_metadata_in)[1] <- "Sample"
inhouse_metadata_in$Sample <- as.character(inhouse_metadata_in$Sample)

## define functions
filter_out_low_expressed <- function(dds, deseq_obj = TRUE){
  library(DESeq2)
  print(paste("no of genes before filtering...", nrow(dds)))
  if (deseq_obj == TRUE){
    # returns a vector of whether the total count of each gene is >= 10 (True or false)
    keep <- rowSums(counts(dds)) >= 10
    # only keep rows (genes) for which keep is TRUE
    dds <- dds[keep,]
    # at least X samples with a count of 10 or more, where X is 5% of samples
    X <- round(0.05*ncol(dds))
    keep <- rowSums(counts(dds) >= 10) >= X
    dds <- dds[keep,]
  } else {
    # returns a vector of whether the total count of each gene is >= 10 (True or false)
    keep <- rowSums(dds) >= 10
    # only keep rows (genes) for which keep is TRUE
    dds <- dds[keep,]
    # at least X samples with a count of 10 or more, where X is 5% of samples
    X <- round(0.05*ncol(dds))
    keep <- rowSums(dds >= 10) >= X
    dds <- dds[keep,]
  }
  print(paste("no of genes after filtering...", nrow(dds)))
  return(dds)
}

### Read in In-house data

inhouse_metadata_combined <- full_join(inhouse_metadata_in, metadata_original_inhouse)
files_inhouse <- file.path(inhouse_metadata_combined$directory, inhouse_metadata_combined$Sample, "quant.sf")
coldata_inhouse <- data.frame(files = files_inhouse, names=inhouse_metadata_combined$Sample, Study = inhouse_metadata_combined$Study, 
                              Subpopulation = inhouse_metadata_combined$Subpopulation, 
                              Tumor_JuxtaTumor = inhouse_metadata_combined$Tumor_JuxtaTumor,
                              Patient = inhouse_metadata_combined$Patient,
                              stringsAsFactors=FALSE)
# tx2gene generated with generate_tx2gene_table.R or supplied
tx2gene <- read.table(tx2gene_file, header = T)
# read in quant.sf file for each sample
se <- tximeta(coldata_inhouse, skipMeta=TRUE, txOut=FALSE, tx2gene=tx2gene)
dds <- DESeqDataSet(se, design = ~1)
dds <- filter_out_low_expressed(dds)
inhouse_counts <- as.data.frame(assays(dds)$counts)
inhouse_counts$Gene <- str_split_fixed(rownames(inhouse_counts), pattern = "\\.", n = 2)[,1]

### Read in GTEx data

gtex_in <- vroom(gtex_file, skip = 2)
gtex_in <- gtex_in %>% dplyr::select(-c(id, Description))
gtex_in <- as.data.frame(gtex_in)
rownames(gtex_in) <- gtex_in$Name
gtex_in$Name <- NULL
gtex_in <- filter_out_low_expressed(gtex_in, deseq_obj = F)
gtex_in$Gene <- str_split_fixed(rownames(gtex_in), pattern = "\\.", n = 2)[,1]
gtex_in[1:5,1:5]
#dds_gtex <- DESeqDataSetFromMatrix(gtex_in)

### Match by genes

#### first get median of duplicated ENSGs

get_median_duplicated_genenames <- function(df){
  n_occur <- data.frame(table(df$Gene))
  dups <- df[df$Gene %in% n_occur$Var1[n_occur$Freq > 1],]
  not_dups <- df[!(df$Gene %in% n_occur$Var1[n_occur$Freq > 1]),]
  dups_summarise <- dups %>% group_by(Gene) %>% dplyr::summarise((across(where(is.numeric), median)))
  df_not_dup <- rbind.data.frame(not_dups, dups_summarise)
  return(df_not_dup)
}
inhouse_counts_not_dup <- get_median_duplicated_genenames(inhouse_counts)
gtex_in_not_dup <- get_median_duplicated_genenames(gtex_in)


#### match by genes

inhouse_gtex_join <- full_join(inhouse_counts_not_dup, gtex_in_not_dup)
rownames(inhouse_gtex_join) <- inhouse_gtex_join$Gene
inhouse_gtex_join$Gene <- NULL
#inhouse_gtex_join <- lapply(na.omit(inhouse_gtex_join), as.numeric)
inhouse_gtex_join <- mutate_all(na.omit(inhouse_gtex_join), function(x) as.integer(as.character(x)))

#### make coldata
coldata_gtex <- data.frame(names = colnames(gtex_in)[1:ncol(gtex_in)-1], Study = rep("GTEX", ncol(gtex_in)-1), Tumor_JuxtaTumor = rep("normal", ncol(gtex_in)-1))
coldata_inhouse_reduced <- coldata_inhouse %>% dplyr::select(c(names, Study, Tumor_JuxtaTumor))
coldata_inhouse_gtex <- rbind.data.frame(coldata_inhouse_reduced, coldata_gtex)

batchQC(inhouse_gtex_join, batch=coldata_inhouse_gtex$Study, 
        report_file="batchqc_report_gtex_fibroblast.html", report_dir=".", 
        report_option_binary="111111111",
        view_report=FALSE, interactive=FALSE, batchqc_output=TRUE)


