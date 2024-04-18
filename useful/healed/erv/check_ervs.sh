#!/bin/bash
set -e
# get BED file for the enrichment panel
wget https://emea.support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/illumina-prep/exome/hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED 
# get list of ervs
wget https://g-a8b222.dd271.03c0.data.globus.org/pub/databases/genenames/hgnc/tsv/locus_types/endogenous_retrovirus.txt
# extract second column of endogenous_retrovirus.txt excluding the header
cut -f 2 endogenous_retrovirus.txt | tail -n +2 > endogenous_retrovirus_in_hgnc.txt
# count the number of ERVs in the panel
echo "Number of ERVs in the panel:"
grep -f endogenous_retrovirus_in_hgnc.txt hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED | wc -l
# write the ERVs in the panel to a file
grep -f endogenous_retrovirus_in_hgnc.txt hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED > ERVs_in_panel.bed
# count the number of ERVs in  not in the panel
ervs=$(cat endogenous_retrovirus_in_hgnc.txt | wc -l)
ervs_in_panel=$(grep -f endogenous_retrovirus_in_hgnc.txt hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED | wc -l)
echo "Number of ERVs not in panel"
echo $((ervs - ervs_in_panel))
# clean up files
rm hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED* endogenous_retrovirus_in_hgnc.txt endogenous_retrovirus.txt*
