#!/bin/bash
set -e
# get BED file for the enrichment panel if not present
if [ ! -f hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED ]; then
    wget https://emea.support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/illumina-prep/exome/hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED 
fi
#wget https://emea.support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/illumina-prep/exome/hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED 
# extract second column of hgnc_non_coding_rnas.txt excluding the header
cut -f 2 hgnc_non_coding_rnas.txt | tail -n +2 > lincs_in_hgnc.txt
# count number of LINCs in lincs_in_hgnc.txt are present in the BED file
echo "number of LINCs in hgnc_non_coding_rnas.txt that are present in the BED file:"
grep -of lincs_in_hgnc.txt hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED | sort | uniq | wc -l
# count the number of LINCs in  not in the panel
lincs=$(cat lincs_in_hgnc.txt | wc -l)
lincs_in_panel=$(grep -of lincs_in_hgnc.txt hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED | wc -l)
echo "Number of non-coding RNAs not in panel"
echo $((lincs_in_panel - lincs))
# write lincs in hgnc to a file
grep -f lincs_in_hgnc.txt hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED > lincs_in_both.txt
# remove unecessary files
rm hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED* lincs_in_hgnc.txt
