#!/usr/bin/Rscript
library(optparse)
library(stringr)
library(dplyr)

# input file format
# Pt-No fastq

option_list = list(
 make_option(c("-f", "--file"), type = "character", default = NULL,
                help = "csv file listing files of interest", metavar = "character"),
make_option(c("-d", "--dna"), type = "character", default = "WES",
                help = "sequencing type used for DNA, must be one of WGS, WXS, WES", metavar = "character")
 #make_option(c("-o", "--outfile"), type = "character", default = "out.txt",
  #          help = "outfile", metavar = "character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (length(opt) < 1) {
  print_help(opt_parser)
  stop("must specify infile at least", call. = FALSE)
}

file <- opt$file
dna <- opt$dna
if (!(dna %in% c("WES", "WGS", "WXS"))){
    print_help(opt_parser)
    stop("dna must be one of WGS, WXS, WES", call. = FALSE)

}
#outfile <- opt$outfile

filein <- read.csv(file, header = F)

start_normal_rna <- "nr-"
start_abnormal_rna <- "ar-"
start_normal_dna <- "nd-"
start_abnormal_dna <- "ad-"

outnames_dna <- c()
samplenames <- str_split_fixed(filein$V1, pattern = "\\_", n = 3)[,1]
for (i in 1:nrow(filein)){
    # read1
    if (grepl(filein$V1[i], pattern = "R1") == TRUE) {
        # abnormal read1
        if (filein$V3[i] == "a"){
            newname <- paste("ad-", samplenames[i], "_1.fastq.gz", sep = "")
            outnames_dna <- c(outnames_dna, newname)
        # normal read1
        } else if (filein$V3[i] == "n") {
            newname <- paste("nd-", samplenames[i], "_1.fastq.gz", sep = "")
            outnames_dna <- c(outnames_dna, newname)
        }
    }
    # read2
    else if (grepl(filein$V1[i], pattern = "R2") == TRUE) {
        # abnormal read2
          if (filein$V3[i] == "a"){
            newname <- paste("ad-", samplenames[i], "_2.fastq.gz", sep = "")
            outnames_dna <- c(outnames_dna, newname)
        } else if (filein$V3[i] == "n") {
            # normal read2
            newname <- paste("nd-", samplenames[i], "_2.fastq.gz", sep = "")
            outnames_dna <- c(outnames_dna, newname)
        }
    } else {
        print("something wrong with for loop")
    }
}

outnames_rna <- c()
for (i in 1:nrow(filein)){
    # read1
    if (grepl(filein$V1[i], pattern = "R1") == TRUE) {
        # abnormal read1
        if (filein$V3[i] == "a"){
            newname <- paste("ar-", samplenames[i], "_1.fastq.gz", sep = "")
            outnames_rna <- c(outnames_rna, newname)
        # normal read1
        } else if (filein$V3[i] == "n") {
            newname <- paste("nr-", samplenames[i], "_1.fastq.gz", sep = "")
            outnames_rna <- c(outnames_rna, newname)
        }
    }
    # read2
    else if (grepl(filein$V1[i], pattern = "R2") == TRUE) {
        # abnormal read2
          if (filein$V3[i] == "a"){
            newname <- paste("ar-", samplenames[i], "_2.fastq.gz", sep = "")
            outnames_rna <- c(outnames_rna, newname)
        } else if (filein$V3[i] == "n") {
            # normal read2
            newname <- paste("nr-", samplenames[i], "_2.fastq.gz", sep = "")
            outnames_rna <- c(outnames_rna, newname)
        }
    } else {
        print("something wrong with for loop")
    }
}


#read1 <- filein[grepl(filein$V1, pattern = "R1") == TRUE,]
#read1_a <- read1[read1$V3 == "a",]
#read1_n <- read1[read1$V3 == "n",]
#read2 <- filein[grepl(filein$V1, pattern = "R2") == TRUE,]
#read2_a <- read2[read2$V3 == "a",]
#read2_n <- read2[read2$V3 == "n",]

# create filenames for rna
#normal_rna_r1 <- paste(start_normal_rna, str_split_fixed(read1_n$V1, pattern = "\\_", n = 3)[,1], "_1.fastq.gz", sep = "")
#abnormal_rna_r1 <- paste(start_abnormal_rna, str_split_fixed(read1_a$V1, pattern = "\\_", n = 3)[,1], "_1.fastq.gz", sep = "")
#normal_rna_r2 <- paste(start_normal_rna, str_split_fixed(read2_n$V1, pattern = "\\_", n = 3)[,1], "_2.fastq.gz", sep = "")
#abnormal_rna_r2 <- paste(start_abnormal_rna, str_split_fixed(read2_n$V1, pattern = "\\_", n = 3)[,1], "_2.fastq.gz", sep = "")
#
# create filenames for dna
#normal_dna_r1 <- paste(start_normal_dna, str_split_fixed(read1_n$V1, pattern = "\\_", n = 3)[,1], "_1.fastq.gz", sep = "")
#abnormal_dna_r1 <- paste(start_abnormal_dna, str_split_fixed(read1_a$V1, pattern = "\\_", n = 3)[,1], "_1.fastq.gz", sep = "")
#normal_dna_r2 <- paste(start_normal_dna, str_split_fixed(read2_n$V1, pattern = "\\_", n = 3)[,1], "_2.fastq.gz", sep = "")
#abnormal_dna_r2 <- paste(start_abnormal_dna, str_split_fixed(read2_a$V1, pattern = "\\_", n = 3)[,1], "_2.fastq.gz", sep = "")

# put together
#dna_out_reads_original <- c(read1$V1, read2$V1)
#dna_out_newnames <- c(normal_dna_r1, normal_dna_r2, abnormal_dna_r1, abnormal_dna_r2)
#samples_order_dna <- str_split_fixed(dna_out_reads_original, pattern = "\\_", n = 3)[,1]
#dna_out_newnames_order <- dna_out_newnames[names(sapply(dna_out_newnames, function(x) { grep(x, samples_order_dna) }))]
#dna_out <- cbind(dna_out_reads_original, dna_out_newnames_order)

#rna_out_reads_original <- c(read1$V1, read2$V1)
#rna_out_newnames <- c(normal_rna_r1, normal_rna_r2, abnormal_rna_r1, abnormal_rna_r2)
#samples_order_rna <- str_split_fixed(rna_out_reads_original, pattern = "\\_", n = 3)[,1]
#rna_out_newnames_order <- names(sapply(rna_out_newnames, function(x) { grep(x, samples_order_rna) }))
#rna_out <- cbind(rna_out_reads_original, rna_out_newnames_order)

dna_out <- cbind(filein$V1, outnames_dna)
colnames(dna_out) <- NULL
rna_out <- cbind(filein$V1, outnames_rna)
colnames(rna_out) <- NULL

filein_headers <- filein
colnames(filein_headers) <- c("file_original", "patient", "abnormal_normal")
samplesheet_dna <- filein_headers
samplesheet_dna$file_lens <- outnames_dna
samplesheet_dna$Sequencing_Method <- dna

samplesheet_rna <- filein_headers
samplesheet_rna$file_lens <- outnames_rna
samplesheet_rna$Sequencing_Method <- "RNA"

samplesheet_out <- rbind.data.frame(samplesheet_dna, samplesheet_rna)
#write.table(rna_out, file = "rna_out.txt", sep = "\t",
#            col.names = F, row.names = F, quote = F)
#write.table(dna_out, file = "dna_out.txt", sep = "\t",
#            col.names = F, row.names = F, quote = F)

write.csv(samplesheet_out, file = "samplesheet.csv", sep = ",",
            col.names = F, row.names = F, quote = F)
