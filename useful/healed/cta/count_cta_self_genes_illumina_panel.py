#!/usr/bin/python3

# written with ChatGPT
# Open and read lines from both files
with open('cta_and_self_antigen.homo_sapiens.gene_list', 'r') as f1, open('hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED', 'r') as f2:
    lines_file1 = set(f1.read().splitlines())
    lines_file2 = f2.read()

# Count lines from file1 that appear in file2
count = sum(1 for line in lines_file1 if line in lines_file2)

# for sanity gather up genes that are not in file2 and print
not_in_file2 = []

for line in lines_file1:
    if not line in lines_file2:
        not_in_file2.append(line)

# percent of CTA/self antigen genes covered by Twist Exome Panel
percent = 100*(count/len(lines_file1))
print("Number of CTA/Self genes found in Illumina Exome 2.5:", count)
print("Percentage of CTA/Self genes found in Illumina Exome 2.5: ", percent)
print("CTA/Self genes not found in Illumina Exome 2.5:", not_in_file2)