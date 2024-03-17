#!/usr/bin/Rscript

# The purpose of this script is to do the analysis for chr2_9986202_G_A_b38 and see if its targets are enriched for trans-eQTLs compared to all genes
# Load in packages
library(stringr)
library(biomaRt)
library(vroom)
library(dplyr, warn.conflicts = FALSE)
library(GenomicRanges, quietly = T)

# Take user input (average expression or use all expression values)
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("At least one argument must be supplied ('average' or 'default' )", call.=FALSE)
} else if (args[1] != "average" && args[1] != "default") {
  stop("Argument must either be 'average' or 'default'", call. = FALSE)
}

print(paste("Arg 1:",args[1], "arg 2:",args[2],"arg 3:", args[3], "arg 4:", args[4], "arg 5:", args[5], sep = " "))
print(paste("Start time:", Sys.time(), sep = ""))

# Point to file names
genotypes <- "/data/kryan/project/gtex/genotype_data_dbdp.vcf"
allele.frequency.desired <- 0.01
barrera.variants <- "/data/kryan/project/gtex/analysis/dbdp_nssnp_barrera_with_GTEx_varid.txt"
if (args[2] == "whole_blood_tpm"){
	expression_file <- "/data/kryan/project/gtex/sample_data/GTEx_Analysis_2017_Whole_blood_tpm.gct"
	cov_file <- "/data/kryan/project/gtex/sample_data/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt"
} else if (args[2] == "whole_blood_norm") {
	expression_file <- "/data/kryan/project/gtex/sample_data/Whole_blood_v8.normalized_filtered_expression.bed"
	cov_file <- "/data/kryan/project/gtex/sample_data/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt"
} else if (args[2] == "frac_rank") {
	#expression_file <- "/data/kryan/project/gtex/analysis/frac_rank.out.all"
	expression_file <- "/data/kryan/project/gtex/sample_data/frac_rank_with_extra_info_v2.txt"
	cov_file <- "/data/kryan/project/gtex/sample_data/covariates_frac_rank_extracted_exp_and_geno.txt"
}
psam <- "/data/kryan/project/gtex/gtex.psam"
dorothea_file <- "/data/kryan/project/gtex/analysis/dorothea_hs_with_target_info.txt"
# Read in genotype table
genotype.data <- read.table(genotypes, stringsAsFactors=FALSE)
# Change genotype encoding from 0/0,0/1,1/1 to 0,1,2
genotype_to_numeric <- function(par){
 par.out <- par
 par.out[which(par == "0/0")] <- 0
 par.out[which(par == "0/1")] <- 1
 par.out[which(par == "1/1")] <- 2
return(as.numeric(par.out))
}

# Read in covariates
if (args[2] == "frac_rank"){
	cov <- read.table(cov_file, check.names=F)
	#cov <- cov[-1,]
	cov_t <- t(cov)

} else if (args[2] == "whole_blood_norm" | args[2] == "whole_blood_tpm") {
	cov <- vroom(cov_file)
	cov_t <- t(cov[-1])
	colnames(cov_t) <- cov$ID
}

if (args[4] == "reduced_covariates"){
        reduced_covariates <- c("PC1", "PC2", "PC3", "PC4", "PC5", "pcr", "platform", "sex")
        cov_t <- cov_t[,reduced_covariates]
}

# convert to numeric - get some warnings
variants.genotypes <- genotype.data[,10:ncol(genotype.data)]
variants.numeric <- apply(variants.genotypes, 2, genotype_to_numeric)
variant.ids <- genotype.data$V3
# Read in expression data
expression_data <- vroom(expression_file)
# Convert sample ids to individual ids if required ("convert_ids" as argument 3)
if (args[3] == "convert_ids"){
	split.up <- str_split_fixed(colnames(expression_data)[3:ncol(expression_data)], "-",5)
	joined.up <- paste(split.up[,1],"-",split.up[,2], sep = "")
	first_2_cols <- colnames(expression_data[1:2])
	all.cols <- c(first_2_cols,joined.up)
	colnames(expression_data) <- all.cols
}

# remove rows with duplicated gene names - we don't know the ensembl ids for the frac_rank expression data so I wanted to do the same to the whole blood data so that we can compare like with like
if (args[2] == "whole_blood_tpm"){
	singletons <- names(which(table(expression_data$Description) == 1))
	expression_data_notduplicated <- expression_data[expression_data$Description %in% singletons, ]
} else if (args[2] == "frac_rank"){
	singletons <- names(which(table(expression_data$Gene) == 1))
	expression_data_notduplicated <- expression_data[expression_data$Gene %in% singletons, ]

} else if (args[2] == "whole_blood_norm"){
	singletons <- names(which(table(expression_data$ensembl_gene_id_target_short) == 1))
	expression_data_notduplicated <- expression_data[expression_data$ensembl_gene_id_target_short %in% singletons, ]
}
# Read in sample ids and match up with the sample ids in the expression data file
psam.ids <- read.table(psam)
colnames(variants.numeric) <- psam.ids$V1
if (args[2] == "whole_blood_tpm"){
	isect <- intersect(colnames(variants.numeric),colnames(expression_data[3,]))
} else if (args[2] == "frac_rank"){
	isect <- intersect(colnames(variants.numeric),colnames(expression_data[7,]))
} else if (args[2] == "whole_blood_norm"){
	isect <- intersect(colnames(variants.numeric),colnames(expression_data[6,]))
}
variants.genotypes.matched <- variants.numeric[,isect]

# Calculate MAFs based on individuals present in data. Exclude NAs from the sum as well as the denominator of the fraction
maf <- function(par){
 par.out <- sum(par, na.rm = TRUE)/(2*length(par[!is.na(par)]))
 par.out <- pmin(par.out, 1-par.out)
 return(par.out)
}
mafs <- apply(variants.genotypes.matched, 1, maf)
#for (k in c(0.05,0.03,0.01)){
#k <- as.numeric(k)
#print(paste("k:", k, sep = ""))
variant.ids.range <- variant.ids[mafs >= allele.frequency.desired]
#variant.ids.range <- variant.ids[mafs >= k]
print(paste("Number of variants present at desired allele frequency:", length(variant.ids.range)))
variants.range.genotypes.matched <- variants.genotypes.matched[which(mafs >= allele.frequency.desired),]
#variants.range.genotypes.matched <- variants.genotypes.matched[which(mafs >= k),]
print(paste("number of rows in variants.range.genotypes.matched", nrow(variants.range.genotypes.matched)))
expression.data.samples.matched <- expression_data_notduplicated[,isect]
if (args[2] == "whole_blood_tpm"){
	rownames(expression.data.samples.matched) <- expression_data_notduplicated$Description
} else if ((args[2] == "frac_rank")){
	rownames(expression.data.samples.matched) <- expression_data_notduplicated$Gene
} else if (args[2] == "whole_blood_norm") {
	rownames(expression.data.samples.matched) <- expression_data_notduplicated$ensembl_gene_id_target_short
}

# Read in table from Barrera paper (with GTEx ids added in manually)
tf.dbdp <- read.table(barrera.variants)

# Only keep rows that are in GTEx at the desired allele frequency
tf.dbdp.variants <- tf.dbdp[tf.dbdp$gtex_var_format_tfdbdp_b38 %in% variant.ids.range,]

# Get rid of duplicated variants (the Barrera paper has separate entries for different transcript ids)
tf.dbdp.uniq <- tf.dbdp.variants[!duplicated(tf.dbdp.variants$gtex_var_format_tfdbdp_b38),]

# Read in dorothea file - created manually with added features, as of 17/06/2021
dorothea <- read.table(dorothea_file, header = T, sep = "\t", check.names = F)

# Only keep the variants in transcription factors that have data in the dorothea database
tf.dbdp.uniq.indorothea <- tf.dbdp.uniq[tf.dbdp.uniq$Gene.symbol %in% dorothea$tf,]

print(paste("Number of transcription factors with GTEx variants in desired frequency range available in dorothea database:",nrow(tf.dbdp.uniq.indorothea)))

# Function for extracting p-value from linear model - don't actually need this anymore
lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}

# For each transcription factor variant, get the targets for the transcription factor and fit a linear model for the expression of the targets wrt increasing numbers of the minor allele
outputs <- list()
no.of.tests <- 0
#df.out <- data.frame(variant = character(), tf.interest = character(), target = character(), p_value = as.numeric())
df.out <- data.frame()
# Create GRanges object from dorothea
df <- data.frame(chrom=dorothea$target_chromosome, start=dorothea$target_start, end=dorothea$target_end)
df$chrom <- paste("chr",df$chrom, sep = "")
gr1 <- as(df, "GRanges")

# Dataframes for dorothea upregulated and downregulated targets
dorothea_upreg <- dorothea[dorothea$mor == 1,]
dorothea_downreg <- dorothea[dorothea$mor == -1,]

# Create GRanges object for filtered variants range
vars <- tf.dbdp.uniq.indorothea[,22]
vars.split <- str_split_fixed(vars, "_",5)
vars.chromosome <- vars.split[,1]
var.positions <- as.numeric(vars.split[,2])
range.start <- (var.positions - 1e06)
range.start[range.start < 1] <- 1
range.end <- as.numeric(var.positions + 1e06)
df2 <- data.frame(chrom = vars.chromosome,start=range.start, end=range.end)
gr2 <- as(df2, "GRanges")

# Create granges object for expression data
if (args[2] == "whole_blood_norm"){
	df3_norm <- data.frame(chrom = expression_data_notduplicated$`#chr`, start = expression_data_notduplicated$start, end = expression_data_notduplicated$end)
	gr3 <- as(df3_norm, "GRanges")
} else if (args[2] == "frac_rank"){
	df3_frac_rank <- data.frame(chrom = expression_data_notduplicated$chromosome, start = expression_data_notduplicated$start, end = expression_data_notduplicated$end)
	df3_frac_rank$chrom <- paste("chr", df3_frac_rank$chrom, sep = "")
	gr3 <- as(df3_frac_rank, "GRanges")
}

# For each transcription factor variant, get the targets for the transcription factor and fit a linear model for the expression of the targets wrt increasing numbers of the minor allele
outputs <- list()
no.of.tests <- 0
#df.out <- data.frame(variant = character(), tf.interest = character(), target = character(), p_value = as.numeric())

outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}
fisher.exact.pvalues <- c()
chisq.pvalues <- c()
no.of.targets <- c()

# function to carry out eqtl analysis on non-targets
non_targets_eqtl <- function (tab, tab_outersect.rm.na, cov_t, p.values.outersect){
		linear.model.outersect <- lm(tab ~ tab_outersect.rm.na$tf.genotypes + as.matrix(cov_t), na.action = na.omit)
                p.values.outersect <- c(p.values.outersect, summary(linear.model.outersect)$coefficients[2,4])
		return(p.values.outersect)
		}

if(args[1] == "default"){
for (i in 1:nrow(tf.dbdp.uniq.indorothea)){
	var <- tf.dbdp.uniq.indorothea[i,22]
	print(paste("Variant number", i, ":", var, sep = " "))
	tf.interest <- tf.dbdp.uniq.indorothea[i,9]	
	tf.genotypes <- variants.range.genotypes.matched[which(variant.ids.range == var),]
	# Remove targets within 1MB from consideration
	dorothea.possible <- dorothea[countOverlaps(gr1,gr2[i]) == 0,]
	# Remove all genes with 1MB from expression file, only looking for transeQTLs
	expression.data.samples.matched.possible <- expression.data.samples.matched[countOverlaps(gr3,gr2[i]) == 0,]
	rownames(expression.data.samples.matched.possible) <- rownames(expression.data.samples.matched)[which((countOverlaps(gr3,gr2[i]) == 0))]
	# Filter so that variants are not on the same chromosome within 1MB of the target gene
	if (args[2] == "whole_blood_norm"){
		tf.targets.interest <- dorothea.possible$ensembl_gene_id_target[dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest]
	} else if (args[2] == "frac_rank" | args[2] == "whole_blood_tpm"){
		tf.targets.interest <- dorothea.possible$target[dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest]
	}
        tf.targets.names <- dorothea.possible$target[dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest]
	mor.targets <- dorothea.possible$mor[which(dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest)]
	confidence.targets <- dorothea.possible$confidence[which(dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest)]
	chrom.targets <- dorothea.possible$target_chromosome[which(dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest)]
	targets.start <- dorothea.possible$target_start[which(dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest)]
	targets.end <- dorothea.possible$target_end[which(dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest)]
	target.range <-paste(targets.start,targets.end, sep = "-" )
	chrom.target.range <- paste(chrom.targets,target.range, sep = ":")
	chrom.target.range <- paste("chr",chrom.target.range, sep = "")
	df.targets <- cbind.data.frame(tf.targets.interest, tf.targets.names, chrom.target.range, mor.targets, confidence.targets)
	df.targets.intersect <- df.targets[which(df.targets$tf.targets.interest %in% rownames(expression.data.samples.matched.possible)),]
	p.values <- c()
	p.values.outersect <- c()
	intersect.tf.targets.expression.data <- intersect(df.targets$tf.targets.interest, rownames(expression.data.samples.matched.possible))
	# outersect introduces some targets that are not possible it seems, introduces as na's in next steps, need to be removed
	outersect.tf.targets.expression.data <- outersect(df.targets$tf.targets.interest, rownames(expression.data.samples.matched.possible))
	expression_data_transposed <- t(expression.data.samples.matched.possible[intersect.tf.targets.expression.data,])
	expression_data_transposed_outersect <- t(expression.data.samples.matched.possible[outersect.tf.targets.expression.data,])
	variants.range.genotypes.matched.vector <- unlist(variants.range.genotypes.matched)
	tab <- cbind.data.frame(tf.genotypes, expression_data_transposed, stringsAsFactors = FALSE)
	tab[,2:ncol(tab)] <- sapply(tab[,2:ncol(tab)], as.numeric)
	tab_outersect <- cbind.data.frame(tf.genotypes, expression_data_transposed_outersect, stringsAsFactors = FALSE)
        tab_outersect[,2:ncol(tab_outersect)] <- sapply(tab_outersect[,2:ncol(tab_outersect)], as.numeric)
	for (j in 2:ncol(tab)){
		if (args[4] == "all_covariates" | args[4] == "reduced_covariates"){
                                linear.model <- lm(tab[,j] ~ tab$tf.genotypes + as.matrix(cov_t), na.action = na.omit)
                                p.values <- c(p.values, summary(linear.model)$coefficients[2,4])
                 } else if (args[4] == "no_covariates"){
                                linear.model <- lm(tab[,j] ~ tab$tf.genotypes, na.action = na.omit)
                                p.values <- c(p.values, summary(linear.model)$coefficients[2,4])
                 }

 		# linear.model <- lm(tab[,j] ~ tab$tf.genotypes + as.matrix(cov_t), na.action = na.omit)
		# p.values <- c(p.values, summary(linear.model)$coefficients[2,4])
		}
	not_all_na <- function(x) any(!is.na(x))
	tab_outersect_rm_na <- tab_outersect %>% select(where(not_all_na))
	#tab_outersect_rm_na_expression <- tab_outersect_rm_na[,2:ncol(tab_outersect)]
	#start.time.eqtl.nontarget <- Sys.time()
	for (l in 2:ncol(tab_outersect_rm_na)){
		#print(l)
		if (args[4] == "all_covariates" | args[4] == "reduced_covariates"){
                                linear.model.outersect <- lm(tab_outersect_rm_na[,l] ~ tab_outersect_rm_na$tf.genotypes + as.matrix(cov_t), na.action = na.omit)
				p.values.outersect <- c(p.values.outersect, summary(linear.model.outersect)$coefficients[2,4])
                 } else if (args[4] == "no_covariates"){
                                linear.model.outersect <- lm(tab_outersect_rm_na[,j] ~ tab_outersect_rm_na$tf.genotypes, na.action = na.omit)
				p.values.outersect <- c(p.values.outersect, summary(linear.model.outersect)$coefficients[2,4])
                 }

		#linear.model.outersect <- lm(tab_outersect_rm_na[,l] ~ tab_outersect_rm_na$tf.genotypes + as.matrix(cov_t), na.action = na.omit)
               #p.values.outersect <- c(p.values.outersect, summary(linear.model.outersect)$coefficients[2,4])
	}
	#p.values.outersect <- sapply(X = tab_outersect_rm_na_expression, non_targets_eqtl, tab_outersect_rm_na = tab_outersect_rm_na, cov_t = cov_t, p.values.outersect = p.values.outersect)
	#end.time.eqtl.nontarget <- Sys.time()
	#print(paste("time taken eqtl nontargets:", end.time.eqtl.nontarget - start.time.eqtl.nontarget))
	eqtl_targets <- length(p.values[p.values < 0.05])
   	eqtl_nontargets <- length(p.values.outersect[p.values.outersect < 0.05])
   	noneqtl_targets <- length(p.values[p.values >= 0.05])
   	noneqtl_nontargets <- length(p.values.outersect[p.values.outersect >= 0.05])
   	contingency_tab <- array(c(eqtl_targets, eqtl_nontargets, noneqtl_targets, noneqtl_nontargets), c(2,2))
	#start.time.fisher <- Sys.time()
	fisher.exact.pvalues <- c(fisher.exact.pvalues, fisher.test(contingency_tab, alternative="greater")$p.value)
	#end.time.fisher <- Sys.time()
	#print(paste("time taken fisher:", end.time.fisher - start.time.fisher))
	chisq.pvalues <- c(chisq.pvalues, chisq.test(contingency_tab)$p.value)
	#results <- cbind(rep(var, length(intersect.tf.targets.expression.data)), 
	 #                rep(tf.interest, length(intersect.tf.targets.expression.data)),
	  #               intersect.tf.targets.expression.data, df.targets.intersect$tf.targets.names, df.targets.intersect$chrom.target.range, p.values, cof, df.targets.intersect$mor.targets, df.targets.intersect$confidence.targets)
	#df.out <- rbind.data.frame(df.out,results)
	no.of.targets <- c(no.of.targets,(ncol(tab)-1))
  }
  FDR.adjusted.pvals.fisher <- p.adjust(fisher.exact.pvalues, method = "BH")
  FDR.adjusted.pvals.chisq <- p.adjust(chisq.pvalues, method = "BH")
  results <- cbind.data.frame(tf.dbdp.uniq.indorothea[,22],tf.dbdp.uniq.indorothea[,9], fisher.exact.pvalues, FDR.adjusted.pvals.fisher, chisq.pvalues, FDR.adjusted.pvals.chisq, no.of.targets)
  colnames(results) <- c("Variant", "Variant_gene", "fisher.exact.pvalues", "FDR.adjusted.pvals.fisher", "chisq.pvalues", "FDR.adjusted.pvals.chisq", "Number_of_targets")
  results.sorted <- results[order(results$FDR.adjusted.pvals.fisher),]
  write.table(results.sorted, file = paste("/data/kryan/project/gtex/outputs/trans_eqtl_results/", args[2], "/all_targets/", args[4], "/", args[5], "/enrichment_", args[2], "_", args[4], "_AF1.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
  print(paste("End time:", Sys.time(), sep = ""))
  # print(paste("Number of genes:", nrow(expression.data.samples.matched.possible)))
  # print(paste("Number of tests - targets:", length(p.values)))
  # print(paste("Number of tests - non-targets:", length(p.values.outersect)))
  # print(paste("total number of tests (should be equal to number of genes):", length(p.values)+length(p.values.outersect)))
  # print(paste("Number of nominally significant trans-eQTLs, targets", length(p.values[p.values < 0.05])))
  # print(paste("Number of non-significant trans-eQTLs, targets", length(p.values[p.values >= 0.05])))
  # print(paste("Number of nominally significant trans-eQTLs, non-targets", length(p.values.outersect[p.values.outersect <= 0.05])))
  # print(paste("Number of non-significant trans-eQTLs, non-targets", length(p.values.outersect[p.values.outersect > 0.05])))

} else {
	df.out.upreg <- data.frame()
	df.out.downreg <- data.frame()
	no.of.tests <- 0
	  for (i in 1:nrow(tf.dbdp.uniq.indorothea)){
	    var <- tf.dbdp.uniq.indorothea[i,22]
	    tf.interest <- tf.dbdp.uniq.indorothea[i,9]	
	    tf.genotypes <- variants.range.genotypes.matched[which(variant.ids.range == var),]
	    # Remove targets within 1MB from consideration
	    dorothea.possible <- dorothea[countOverlaps(gr1,gr2[i]) == 0,]
	    # Filter so that variants are not on the same chromosome within 1MB of the target gene
	    if (args[2] == "whole_blood_norm"){
                tf.targets.interest <- dorothea.possible$ensembl_gene_id_target[dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest]
            } else if (args[2] == "frac_rank" | args[2] == "whole_blood_tpm"){
                tf.targets.interest <- dorothea.possible$target[dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest]
            }
	    mor.targets <- dorothea$mor[which(dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest)]
	    confidence.targets <- dorothea$confidence[which(dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest)]
	    chrom.targets <- dorothea.possible$target_chromosome[which(dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest)]
            targets.start <- dorothea.possible$target_start[which(dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest)]
            targets.end <- dorothea.possible$target_end[which(dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest)]
            target.range <-paste(targets.start,targets.end, sep = "-" )
            chrom.target.range <- paste(chrom.targets,target.range, sep = ":")
            chrom.target.range <- paste("chr",chrom.target.range, sep = "")
	    tf.targets.upreg <- tf.targets.interest[mor.targets == 1]
	    tf.targets.downreg <- tf.targets.interest[mor.targets == -1]
	    confidence.targets.upreg <- confidence.targets[mor.targets == 1]
	    confidence.targets.downreg <- confidence.targets[mor.targets == -1]
	    p.values.upreg <- c()
	    p.values.downreg <- c()
	    cof.upreg <- c()
	    cof.downreg <- c()
	    variants.range.genotypes.matched.vector <- unlist(variants.range.genotypes.matched)
	    if (length(tf.targets.upreg) == 0){
	    	cof.upreg <- NA
            	p.values.upreg <- NA
		print(paste("No upregulated genes in variant:", i))

            } else {
	    	intersect.tf.targets.upreg.expression.data <- intersect(tf.targets.upreg, rownames(expression.data.samples.matched))
            	expression_data_transposed_upreg <- t(expression.data.samples.matched[intersect.tf.targets.upreg.expression.data,])
            	tab_upreg <- cbind.data.frame(tf.genotypes, expression_data_transposed_upreg, stringsAsFactors = FALSE)
		#print(paste("dim(tab_upreg)", i, dim(tab_upreg)))
            	if (ncol(tab_upreg) == 2){
            #tab_rank_upreg <- apply(tab_upreg[,2], FUN = rank, MARGIN = 2)
            #tab_rank_upreg <- cbind(tab_upreg$tf.genotypes, tab_rank_upreg)
            #tab_rank_upreg_mean <- tab_rank_upreg$expression_data_transposed_upreg
	    #names(tab_rank_upreg_mean) <- rownames(tib_rank_upreg)
	    		tab_rank_upreg <- rank(tab_upreg[,2])
            #tab_rank_downreg <- cbind(tab_downreg$tf.genotypes, tab_rank_downreg)
            #tab_rank_downreg_mean <- tab_rank_downreg$expression_data_transposed_downreg
            		tab_rank_upreg_mean <- tab_rank_upreg

                } else {
	    		tab_rank_upreg <- apply(tab_upreg[,2:ncol(tab_upreg)], FUN = rank, MARGIN = 2)
	    		tab_rank_upreg <- cbind(tab_upreg$tf.genotypes, tab_rank_upreg)
	    		tab_rank_upreg_mean <- apply(tab_rank_upreg[,2:ncol(tab_rank_upreg)], FUN = mean, MARGIN = 1)
		}
   	    	tab_rank_upreg_mean_genotypes <- cbind.data.frame(tab_upreg$tf.genotypes, tab_rank_upreg_mean)
            	colnames(tab_rank_upreg_mean_genotypes)[1] <- "tf.genotypes"
            	linear_model_upreg <- lm(tab_rank_upreg_mean_genotypes$tab_rank_upreg_mean ~ tab_rank_upreg_mean_genotypes$tf.genotypes + as.matrix(cov_t), na.action = na.exclude)
            	cof.upreg <- c(cof.upreg,summary(linear_model_upreg)$coefficients[2,1])
            	p.values.upreg <- c(p.values.upreg, summary(linear_model_upreg)$coefficients[2,4])
	    	no.of.tests <- no.of.tests + 1
            	print(paste("Upreg,", i))
		}
            
            if (length(tf.targets.downreg) == 0){
	    	cof.downreg <- NA
            	p.values.downreg <- NA
            print(paste("No downregulated target genes for variant:", i))
	    } else {
            	print(paste("Trying downreg:", i))
	    	intersect.tf.targets.downreg.expression.data <- intersect(tf.targets.downreg, rownames(expression.data.samples.matched))
		if (length(intersect.tf.targets.downreg.expression.data) == 0){
		print("No overlap with expression data")
		cof.downreg <- NA
                p.values.downreg <- NA
		} else {
            	expression_data_transposed_downreg <- t(expression.data.samples.matched[intersect.tf.targets.downreg.expression.data,])
            	tab_downreg <- cbind.data.frame(tf.genotypes, expression_data_transposed_downreg, stringsAsFactors = FALSE)
	    	if (ncol(tab_downreg) == 2){
	    		tab_rank_downreg <- rank(tab_downreg[,2])
	   		tab_rank_downreg_mean <- tab_rank_downreg
		} else {
	    		tab_rank_downreg <- apply(tab_downreg[,2:ncol(tab_downreg)], FUN = rank, MARGIN = 2)
	    		tab_rank_downreg <- cbind(tab_downreg$tf.genotypes, tab_rank_downreg)
            		tab_rank_downreg_mean <- apply(tab_rank_downreg[,2:ncol(tab_rank_downreg)], FUN = mean, MARGIN = 1)
	    	}
            	tab_rank_downreg_mean_genotypes <- cbind.data.frame(tab_downreg$tf.genotypes, tab_rank_downreg_mean)
            	colnames(tab_rank_downreg_mean_genotypes)[1] <- "tf.genotypes"
            	linear_model_downreg <- lm(tab_rank_downreg_mean_genotypes$tab_rank_downreg_mean ~ tab_rank_downreg_mean_genotypes$tf.genotypes + as.matrix(cov_t), na.action = na.exclude)
            	cof.downreg <- c(cof.downreg,summary(linear_model_downreg)$coefficients[2,1])
            	p.values.downreg <- c(p.values.downreg, summary(linear_model_downreg)$coefficients[2,4])
            	no.of.tests <- no.of.tests + 1
	    	print(paste("Downreg,", i))
		}
            }
		  #tab_rank_mean_genotypes <- cbind.data.frame(tab$tf.genotypes, tab_rank_mean)
		  #colnames(tab_rank_mean_genotypes)[1] <- "tf.genotypes" 
		  #linear.model <- lm(tab_rank_mean_genotypes$tab_rank_mean ~ tab_rank_mean_genotypes$tf.genotypes)
		 # p.values <- c(p.values, lmp(linear.model))
	    results.upreg <- cbind(var, tf.interest, p.values.upreg, cof.upreg)
	    results.downreg <- cbind(var,tf.interest, p.values.downreg, cof.downreg)
       	    df.out.upreg <- rbind.data.frame(df.out.upreg,results.upreg)
	    df.out.downreg <- rbind.data.frame(df.out.downreg, results.downreg)
	}
	#colnames(df.out) <-  c("Variant_GTEx_ID", "TF_containing_variant", "p_value_upreg", 
		                     #  "p_value_downreg", "coefficient_upreg", "coefficient_downreg")
	colnames(df.out.upreg) <-  c("Variant_GTEx_ID", "TF_containing_variant", "p_value_upreg", "coefficient_upreg")
	colnames(df.out.downreg) <- c("Variant_GTEx_ID", "TF_containing_variant", "p_value_downreg", "coefficient_downreg")
        print(paste("Number of tests carried out", no.of.tests))
        sig.values.df.out <- 0.05/no.of.tests
        print(paste("Significance value:", sig.values.df.out))
        #colnames(df.out) <- c(colnames(results),"Mode_of_regulation","Target_confidence")
        df.out.upreg$"p_value" <- as.numeric(df.out.upreg$"p_value")
	df.out.downreg$"p_value" <- as.numeric(df.out.downreg$"p_value")
        #df.out.na.remove <- na.omit(df.out)
	df.out.upreg.na.remove <- na.omit(df.out.upreg)
	df.out.downreg.na.remove <- na.omit(df.out.downreg)
        print("Output - average method")
	print(df.out.upreg.na.remove, row.names = FALSE)
	print(df.out.downreg.na.remove, row.names = FALSE)
        #sig.results <- df.out[df.out.na.remove$"p_value" <= sig.values.df.out,]
	sig.results.upreg <- df.out.upreg.na.remove[df.out.upreg.na.remove$"p_value" <= sig.values.df.out,]
	sig.results.downreg <- df.out.downreg.na.remove[df.out.downreg.na.remove$"p_value" <= sig.values.df.out,]
	print("Significant results upreg (Bonferroni)")
	print(sig.results.upreg, row.names = FALSE)
	print("Significant results downreg (Bonferroni)")
	print(sig.results.downreg, row.names = FALSE)
        #all.pvals <- c(df.out.na.remove$"p_value_upreg",df.out.na.remove$"p_value_downreg")
        all.pvals <- c(df.out.upreg.na.remove$"p_value_upreg",df.out.downreg.na.remove$"p_value_downreg")
        print(all.pvals)
        FDR.adjusted.pvals <- p.adjust(as.numeric(all.pvals), method = "BH")
        print("FDR adjusted pvals")
        print(FDR.adjusted.pvals, quote = F)
        #FDR.adjusted.pvals.upreg <- FDR.adjusted.pvals[1:length(!is.na(FDR.adjusted.pvals))/2]
        #FDR.adjusted.pvals.downreg <- FDR.adjusted.pvals[length(!is.na(FDR.adjusted.pvals))/2:length(!is.na(FDR.adjusted.pvals))]
	FDR.adjusted.pvalues.upreg <- FDR.adjusted.pvals[1:length(df.out.upreg.na.remove$"p_value_upreg")]
	FDR.adjusted.pvalues.downreg <- FDR.adjusted.pvals[((length(FDR.adjusted.pvalues.upreg))+1):length(FDR.adjusted.pvals)]
        #FDR.adjusted.results <- cbind(df.out, FDR.adjusted.pvals.upreg, FDR.adjusted.pvals.downreg)
	print(paste("length of FDR.adjusted.pvalues.upreg:", length(FDR.adjusted.pvalues.upreg)))
	print(paste("number of rows in df.out.upreg.na.remove:", nrow(df.out.upreg.na.remove)))
	print(paste("length of FDR.adjusted.pvalues.downreg:", length(FDR.adjusted.pvalues.downreg)))
        print(paste("number of rows in df.out.downreg.na.remove:", nrow(df.out.downreg.na.remove)))
        FDR.adjusted.results.upreg <- cbind(df.out.upreg.na.remove, FDR.adjusted.pvalues.upreg)
	FDR.adjusted.results.downreg <- cbind(df.out.downreg.na.remove, FDR.adjusted.pvalues.downreg)
	print("FDR adjusted results upreg")
        print(FDR.adjusted.results.upreg, row.names = FALSE)
	print("FDR adjusted results downreg")
	print(FDR.adjusted.results.downreg, row.names = FALSE)
	sig.results.bh.upreg <- FDR.adjusted.results.upreg[as.numeric(FDR.adjusted.results.upreg$FDR.adjusted.pvalues.upreg) <= 0.05,]
	sig.results.bh.downreg <- FDR.adjusted.results.downreg[as.numeric(FDR.adjusted.results.downreg$FDR.adjusted.pvalues.downreg) <= 0.05,]
	#sig.results.bh.added.upreg <- sig.results.bh.upreg %>% mutate(FDR.adjusted.results.upreg[FDR.adjusted.results.upreg$FDR.adjusted.pvals.upreg <= 0.05])
	#sig.results.bh.added.downreg <- sig.results.bh.downreg %>% mutate(FDR.adjusted.results.downreg[FDR.adjusted.results.downreg$FDR.adjusted.pvals.downreg <= 0.05])
	print("FDR adjusted significant results upreg") 
	print(sig.results.bh.upreg, row.names = FALSE)
	print("FDR adjusted significant results downreg")
	print(sig.results.bh.downreg, row.names = FALSE)
	#write.table(FDR.adjusted.results.upreg, file = paste("/data/kryan/project/gtex/outputs/trans_eqtl_results/whole_blood_results/avg_targets/22062021/",args[2],"_FDR_adjusted_results_upreg_AF",(j*100),".txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
	#write.table(FDR.adjusted.results.downreg, file = paste("/data/kryan/project/gtex/outputs/trans_eqtl_results/whole_blood_results/avg_targets/22062021/",args[2],"_FDR_adjusted_results_downreg_AF",(j*100),".txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
	#write.table(sig.results.upreg, file = paste("/data/kryan/project/gtex/outputs/trans_eqtl_results/whole_blood_results/avg_targets/22062021/",args[2],"_bonferroni_significant_results_upreg_AF",(j*100),".txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
	#write.table(sig.results.downreg, file = paste("/data/kryan/project/gtex/outputs/trans_eqtl_results/whole_blood_results/avg_targets/22062021/",args[2],"_bonferroni_significant_results_downreg_AF",(j*100),".txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
	#write.table(sig.results.bh.upreg, file = paste("/data/kryan/project/gtex/outputs/trans_eqtl_results/whole_blood_results/avg_targets/22062021/",args[2],"_bh_significant_results_upreg_AF",(j*100),".txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
	#write.table(sig.results.bh.downreg, file = paste("/data/kryan/project/gtex/outputs/trans_eqtl_results/whole_blood_results/avg_targets/22062021/",args[2],"_bh_significant_results_downreg_AF",(j*100),".txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)

	}
#}
# Print the bonferroni threshold
#bonferroni.cutoff <- 0.05/no.of.tests
#print(paste("Number of tests carried out", nrow(df.out)))
#sig.values.df.out <- 0.05/nrow(df.out)
#print(paste("Significance value:", sig.values.df.out))
#colnames(df.out) <- c(colnames(results),"Mode_of_regulation","Target_confidence")
#jpeg("/data/kryan/project/gtex/outputs/trans_eqtl_results/whole_blood_trans_eqtl_analysis_hist_pvals_nocovariates.jpeg")
#hist(as.numeric(df.out$"p-value"), xlab = "Unadjusted p-values", main = "eQTL analysis p-values with Whole Blood GTEx Data")
#dev.off()
#print("results table:")
#print(df.out, row.names = FALSE)
#df.out$"p-value" <- as.numeric(df.out$"p_value")
#df.out.na.remove <- na.omit(df.out)
#sig.results <- df.out.na.remove[df.out.na.remove$"p_value" <= sig.values.df.out,]
#write.table(df.out,"/data/kryan/project/gtex/outputs/trans_eqtl_results/whole_blood_results/whole_blood_transeqtl_analysis_bonf05_no_covariates_10062021_summarised_expression_allresults_AF03.txt",quote = FALSE, sep = "\t", row.names = FALSE)
#print("Significant results (Bonferroni)")
#print(sig.results, row.names = FALSE)
#FDR.adjusted.pvals <- p.adjust(df.out.na.remove$"p-value", method = "BH") # adjusted using Benjamini and Hochberg method

#FDR.adjusted.results <- cbind(df.out.na.remove, FDR.adjusted.pvals)
#print(FDR.adjusted.results, row.names = FALSE)
#sig.results.bh <- df.out.na.remove[FDR.adjusted.pvals <= 0.05,]
#sig.results.bh.added <- sig.results.bh %>%
                      #  mutate(FDR.adjusted.pvalues = FDR.adjusted.pvals[FDR.adjusted.pvals <= 0.05]) print("FDR adjusted results") print(sig.results.bh.added, row.names = FALSE) 
#write.table(FDR.adjusted.results, "/data/kryan/project/gtex/outputs/trans_eqtl_results/whole_blood_results/all_targets/whole_blood_trans_eqtl_analysis_no_covariates_allresults_AF05_11062021.txt", quote = FALSE, sep = "\t", row.names = FALSE) 
#write.table(sig.results, "/data/kryan/project/gtex/outputs/trans_eqtl_results/whole_blood_results/all_targets/whole_blood_trans_eqtl_analysis_no_covariates_sigresults_bonf_AF05_11062021.txt", quote = FALSE, sep = "\t", row.names = FALSE) 
#write.table(sig.results.bh.added, "/data/kryan/project/gtex/outputs/trans_eqtl_results/whole_blood_results/all_targets/whole_blood_trans_eqtl_analysis_no_covariates_sigresults_bh_AF05_11062021.txt", quote = FALSE, sep = "\t", row.names = FALSE)

