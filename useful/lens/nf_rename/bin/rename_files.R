#!/usr/bin/Rscript
library(optparse)
library(stringr)

option_list = list(
 make_option(c("-f", "--file"), type="character", default=NULL, help="txt file listing files of interest", metavar="character"),
 make_option(c("-o", "--outfile"), type="character", default="out.txt", help="outfile", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (length(opt) < 1){
  print_help(opt_parser)
  stop("must specify infile", call.=FALSE)
}

opt$file -> file
opt$outfile -> outfile

filein <- read.table(file)
unique_symbols <- str_split_fixed(filein$V1, n = 3, pattern = "\\.")
samplenames <- str_split_fixed(unique_symbols[,1], n = 4, pattern = "\\_")[,1]
patient <- "06"
run_name_start <- c("ad", "ar", "ad", "ar", "nd", "nr", "nd", "nr")
read_pair <- str_split_fixed(unique_symbols[,1], n = 4, pattern = "\\_")[,2]
read_pair <- substr(read_pair, 2, 2)
newname_start <- paste(run_name_start, "Pt", patient, samplenames, sep = "-")
newname_mid <- paste(patient, samplenames, "_", read_pair, sep = "")
newname_end <- "fastq.gz"
newname_full <- paste(newname_start, newname_mid, newname_end, sep = ".")
newtable <- data.frame(filein$V1, newname_full)
colnames(newtable) <- NULL
write.table(newtable, file = outfile, sep = "\t", quote = F, row.names = F)
