#!/usr/bin/python3
# script to count the number of genes in the TIL10 signature that are in the Twist Exome Panel
import pdb
import pandas as pd

# Read TIL10 signature into a DataFrame
df1 = pd.read_table("TIL10_signature.txt")
# Read file 2 into a DataFrame, extracting the relevant column
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

# Find common gene symbols between file 1 and file 2
common_gene_symbols = set(df1['ID']).intersection(gene_symbols_in_file2)
not_common_gene_symbols = set(df1['ID']).difference(gene_symbols_in_file2)
not_common_gene_symbols2 = set(df1['ID']).difference(common_gene_symbols)
# Print the count of common gene symbols
print("Number of genes in TIL10 signature found in Illumina Twist Exome kit:", len(common_gene_symbols))
print("Number of genes in TIL10 signature not found in Illumina Twist Exome kit:", len(not_common_gene_symbols))
print("Genes in TIL10 signature not found in file 2:", not_common_gene_symbols2)
print("Percentage of genes in TIL10 signature found in Illumina Twist Exome kit:", round(len(common_gene_symbols)/len(df1['ID'])*100, 1))

# Read the LM22 signature into a DataFrame
df3 = pd.read_table("LM22.txt")
# Extract the gene symbols from the LM22 signature
lm22_gene_symbols = set(df3['Gene symbol'])
# Find common gene symbols between LM22 and file 2
common_gene_symbols_lm22 = lm22_gene_symbols.intersection(gene_symbols_in_file2)
not_common_gene_symbols_lm22 = lm22_gene_symbols.difference(gene_symbols_in_file2)
# Print the count of common gene symbols
print("Number of genes in LM22 signature found in Illumina Twist Exome kit:", len(common_gene_symbols_lm22))
print("Number of genes in LM22 signature not found in Illumina Twist Exome kit:", len(not_common_gene_symbols_lm22))
print("Percentage of genes in LM22 signature found in Illumina Twist Exome kit:", round(len(common_gene_symbols_lm22)/len(lm22_gene_symbols)*100, 1))
print( "Genes in LM22 signature not found in file 2:", not_common_gene_symbols_lm22)

print("reading in xcell signatures")
df5 = pd.read_csv("xcell_signatures.csv", sep = ",")
# view head of the dataframe
print(df5.head())
gene_symbols_in_file_5 = set()
# parse each column from third column to the end
# if a gene symbol is not already in gene_symbols_in_file_5 and is not na, then add it to the set
for column in df5.columns[2:]:
    for gene in df5[column]:
        if not str(gene) in gene_symbols_in_file_5 and not pd.isna(gene):
            gene_symbols_in_file_5.add(gene)
# Find common gene symbols between file 5 and file 2
common_gene_symbols_xcell = gene_symbols_in_file_5.intersection(gene_symbols_in_file2)
not_common_gene_symbols_xcell = gene_symbols_in_file_5.difference(gene_symbols_in_file2)
# Print the count of common gene symbols
print("Number of genes in xcell signatures found in Illumina Twist Exome kit:", len(common_gene_symbols_xcell))
print("Number of genes in xcell signatures not found in Illumina Twist Exome kit:", len(not_common_gene_symbols_xcell))
print("Percentage of genes in xcell signatures found in Illumina Twist Exome kit:", round(len(common_gene_symbols_xcell)/len(gene_symbols_in_file_5)*100, 1))
print("Genes in xcell signatures not found in file 2:", not_common_gene_symbols_xcell)

# Read the mcp counter signature into a DataFrame
df6 = pd.read_table("mcpcounter_genes.txt", usecols=[0])
gene_symbols_in_file_6 = set()
# Extract the gene symbols from the mcp counter signature
for gene in df6['HUGO symbols']:
    if not str(gene) in gene_symbols_in_file_6:
        gene_symbols_in_file_6.add(gene)
# Find common gene symbols between mcp counter signature and file 2
common_gene_symbols_mcp = gene_symbols_in_file_6.intersection(gene_symbols_in_file2)
not_common_gene_symbols_mcp = gene_symbols_in_file_6.difference(gene_symbols_in_file2)
# Print the count of common gene symbols
print("Number of genes in mcp counter signature found in Illumina Twist Exome kit:", len(common_gene_symbols_mcp))
print("Number of genes in mcp counter signature not found in Illumina Twist Exome kit:", len(not_common_gene_symbols_mcp))
print("Percentage of genes in mcp counter signature found in Illumina Twist Exome kit:", round(len(common_gene_symbols_mcp)/len(gene_symbols_in_file_6)*100, 1))
print("Genes in mcp counter signature not found in file 2:", not_common_gene_symbols_mcp)
