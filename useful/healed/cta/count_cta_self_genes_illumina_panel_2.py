#!/usr/bin/python3

# Open and read lines from both files

import pandas as pd

# Read CTA/self genes into a DataFrame
df1 = pd.read_table("cta_and_self_antigen.homo_sapiens.gene_list", names = ["ID"])
# Read BED file into a DataFrame, extracting the relevant column
df2 = pd.read_table('hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED', sep='\t', header=None, usecols=[3])
df2 = df2.dropna(subset=[3])

# Flatten the entries in the fourth column of df2 and extract unique gene symbols
gene_symbols_in_file2 = set()
for entry in df2[3]:
    #breakpoint()
    #print("entry:", entry)
    gene_info = entry.split(',')
    for gene in gene_info[:2]:
        if not str(gene) in gene_symbols_in_file2:
          # Consider only the first two comma-delimited fields
            gene_symbols_in_file2.add(gene)

common_gene_symbols = set(df1['ID']).intersection(gene_symbols_in_file2)
not_common_gene_symbols = set(df1['ID']).difference(gene_symbols_in_file2)
not_common_gene_symbols2 = set(df1['ID']).difference(common_gene_symbols)

print("Number of unique CTA/self antigen genes used in LENS:", len(set(df1['ID'])))
print("Number of CTA/self antigen genes found in Illumina Twist Exome kit:", len(common_gene_symbols))
print("Number of CTA/self antigen genes not found in Illumina Twist Exome kit:", len(not_common_gene_symbols))
print("CTA/self antigen genes not found in file 2:", not_common_gene_symbols2)
print("Percentage of CTA/self antigen genes found in Illumina Twist Exome kit:", round(len(common_gene_symbols)/len(set(df1['ID']))*100, 1))