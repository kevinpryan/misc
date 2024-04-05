#!/usr/bin/python3
# script to count the number of genes in the LM22 signature that are in the Twist Exome Panel
import pandas as pd
twist_exome = pd.read_table("hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED", header=None, sep="\t", names=["chr", "start", "end", "gene"])
twist_exome = twist_exome.dropna(subset=["gene"])
genes_twist_exome = list(twist_exome["gene"])
lm22 = pd.read_table("LM22.txt")
lm22_genes = list(lm22["Gene symbol"])  

# need to find elements of lm22_genes that are in any part of each element of genes_twist_exome (i.e. grep would return a match)
# and then count them
# for each element of lm22_genes, check if it is in any element of genes_twist_exome
count = 0
count_not_in = 0
count_lines = 1
number_of_lm22_genes_missing = len(lm22_genes)
genes_missing = list(lm22["Gene symbol"]) 
for line in lm22_genes:
    for gene in genes_twist_exome:
        if line in gene:
            count += 1
            number_of_lm22_genes_missing -= 1
            genes_missing.remove(line)
            break
        count_not_in += 1
print(count)
print("number of genes in lm22:", len(lm22_genes))
print("number of genes not in lm22:", number_of_lm22_genes_missing)
print("percentage of lm22 genes in twist exome kit:", count/len(lm22_genes)*100)
print("genes missing:", genes_missing)
#print(count_in)
