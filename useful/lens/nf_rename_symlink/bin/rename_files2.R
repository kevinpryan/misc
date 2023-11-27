#!/usr/bin/Rscript
library(optparse)
library(stringr)
library(dplyr)

option_list = list(
 make_option(c("-f", "--file"), type = "character", default = NULL,
                help = "txt file listing files of interest", metavar = "character"),
 make_option(c("-o", "--outfile"), type = "character", default = "out.txt",
            help = "outfile", metavar = "character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (length(opt) < 1) {
  print_help(opt_parser)
  stop("must specify infile", call. = FALSE)
}

file <- opt$file
outfile <- opt$outfile

filein <- read.table(file)
print(filein$V1)
start_rna <- filein$V1[grepl("r", filein$V1) == TRUE]  %>%
                gsub(., pattern = "Pt-06-", replacement = "") %>%
                strsplit("\\.") %>% lapply(`[[`, 1) %>%
                unlist() %>%
                paste0(., "rna")

#print("start_rna")
#print(start_rna)
start_dna <- filein$V1[grepl("d", filein$V1) == TRUE] %>% 
            gsub(., pattern = "Pt-06-", replacement = "") %>% 
            strsplit("\\.") %>%
            lapply(`[[`, 1) %>%
            unlist() %>%
            paste0(., "dna")

mid_rna <- filein$V1[grepl("r", filein$V1) == TRUE] %>%
            strsplit("\\.") %>%
            lapply(`[[`, 2) %>%
            unlist() %>%
            strsplit("_") %>%
            lapply(`[[`, 1) %>%
            unlist() %>%
            paste0(., "rna")

mid_dna <- filein$V1[grepl("d", filein$V1) == TRUE] %>%
            strsplit("\\.") %>%
            lapply(`[[`, 2) %>%
            unlist() %>%
            strsplit("_") %>%
            lapply(`[[`, 1) %>%
            unlist() %>%
            paste0(., "dna")

end_rna <- as.character(filein$V1[grepl("r", filein$V1) == TRUE]) %>%
            strsplit("_") %>%
            lapply(`[[`, 2) %>%
            unlist()

end_dna <- as.character(filein$V1[grepl("d", filein$V1) == TRUE]) %>% 
            strsplit("_") %>%
            lapply(`[[`, 2) %>%
            unlist()
left_dna <- paste(start_dna, mid_dna, sep = ".")
full_dna <- paste(left_dna, end_dna, sep = "_")
left_rna <- paste(start_rna, mid_rna, sep = ".")
full_rna <- paste(left_rna, end_rna, sep = "_")
full_names <- c(full_rna, full_dna)

files_sorted <- sort(filein$V1)
full_names_sorted <- sort(full_names)
out.table <- cbind(files_sorted, full_names_sorted)
write.table(out.table, file = "out.txt", sep = "\t",
            col.names = F, row.names = F, quote = F)
