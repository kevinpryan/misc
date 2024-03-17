#!/usr/bin/Rscript
library(optparse)

# input file format
# no header
# fastq,Patient_number_as_int,a/n

option_list = list(
 make_option(c("-f", "--file"), type = "character", default = NULL,
                help = "csv file listing files of interest", metavar = "character"),
make_option(c("-d", "--dna"), type = "logical", default = TRUE,
                help = "is DNA data included"),
make_option(c("-r", "--rna"), type = "logical", default = TRUE,
                help = "is RNA data included"),                
make_option(c("-t", "--type_dna"), type = "character", default = "WES",
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

library(stringr)
library(dplyr)

file <- opt$file
if (opt$dna == TRUE){
    dna <- opt$type_dna
    #opt_dna <- TRUE
    if (!(dna %in% c("WES", "WGS", "WXS"))){
    print_help(opt_parser)
    stop("dna must be one of WGS, WXS, WES", call. = FALSE)
    start_normal_dna <- "nd-"
    start_abnormal_dna <- "ad-"
} else {
    print("opt$dna is false")
    #opt_dna <- FALSE
}
}

if (opt$rna == TRUE){
    #opt_rna <- TRUE
    start_normal_rna <- "nr-"
    start_abnormal_rna <- "ar-"
} #else {
  #  opt_rna <- FALSE
#}

#if (opt_rna == FALSE & opt_dna == FALSE){
#  print_help(opt_parser)
#  stop("must specify at least opt_rna or opt_dna as TRUE", call. = FALSE)
#}
if (opt$rna == FALSE & opt$dna == FALSE){
  print_help(opt_parser)
  stop("must specify at least opt_rna or opt_dna as TRUE", call. = FALSE)
}
#outfile <- opt$outfile

filein <- read.csv(file, header = F)

samplenames <- str_split_fixed(filein$V1, pattern = "\\_", n = 3)[,1]
if (opt$dna == TRUE){
outnames_dna <- c()
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
dna_out <- cbind(filein$V1, outnames_dna)
colnames(dna_out) <- NULL
}

if(opt$rna == TRUE) {
outnames_rna <- c()
for (i in 1:nrow(filein)){
    # read1
    #print("filename...")
    #print(filein$V1[i])
    if (grepl(filein$V1[i], pattern = "R1") == TRUE) {
        # abnormal read1
        if (filein$V3[i] == "a"){
            #print("abnormal read 1")
            newname <- paste("ar-", samplenames[i], "_1.fastq.gz", sep = "")
            outnames_rna <- c(outnames_rna, newname)
        # normal read1
        } else if (filein$V3[i] == "n") {
            newname <- paste("nr-", samplenames[i], "_1.fastq.gz", sep = "")
            outnames_rna <- c(outnames_rna, newname)
        } else {
            print("something wrong with for loop read 1")
        }
    }
    # read2
    else if (grepl(filein$V1[i], pattern = "R2") == TRUE) {
        # abnormal read2
        #print("read 2 detected  ")
        #print("normal or abnormal...")
        print(filein$V3[i])
          if (filein$V3[i] == "a"){
            #print("abnormal read 2")
            #print("samplename...")
            #print(samplenames[i])
            newname <- paste("ar-", samplenames[i], "_2.fastq.gz", sep = "")
            #print("newname...")
            #print(newname)
            outnames_rna <- c(outnames_rna, newname)
        } else if (filein$V3[i] == "n") {
            # normal read2
            newname <- paste("nr-", samplenames[i], "_2.fastq.gz", sep = "")
            outnames_rna <- c(outnames_rna, newname)
        }
        else {
            print("something wrong with for loop read 2")
        }
    } else {
        print("something wrong with for loop")
    }
}
rna_out <- cbind(filein$V1, outnames_rna)
colnames(rna_out) <- NULL
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


filein_headers <- filein
colnames(filein_headers) <- c("file_original", "patient", "abnormal_normal")

if (opt$dna == TRUE) {
samplesheet_dna <- filein_headers
samplesheet_dna$file_lens <- outnames_dna
samplesheet_dna$Sequencing_Method <- dna
}

if (opt$rna == TRUE) {
samplesheet_rna <- filein_headers
samplesheet_rna$file_lens <- outnames_rna
samplesheet_rna$Sequencing_Method <- "RNA"
}

if (opt$rna == TRUE & opt$dna == TRUE){
samplesheet_out <- rbind.data.frame(samplesheet_dna, samplesheet_rna)
} else if (opt$rna == TRUE & opt$dna == FALSE) {
    samplesheet_out <- samplesheet_rna
} else if (opt$rna == FALSE & opt$dna == TRUE) {
    samplesheet_out <- samplesheet_dna
}
#write.table(rna_out, file = "rna_out.txt", sep = "\t",
#            col.names = F, row.names = F, quote = F)
#write.table(dna_out, file = "dna_out.txt", sep = "\t",
#            col.names = F, row.names = F, quote = F)

write.csv(samplesheet_out, file = "samplesheet.csv", sep = ",",
            col.names = F, row.names = F, quote = F)
