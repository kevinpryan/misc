# Comparison of TruSeq RNA exome and Twist exome panel

#!/bin/bash
set -e
# get BED file for the Twist enrichment panel if not present
if [ ! -f hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED ]; then
    wget https://emea.support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/illumina-prep/exome/hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED 
fi
# get BED file of TruSeq RNA exome not present
if [ ! -f truseq-rna-exome-targeted-regions-manifest-v1-2.bed ]; then
    wget https://emea.support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/truseq/truseq-rna-exome-targeted-regions-manifest-v1-2-bed.zip
    unzip truseq-rna-exome-targeted-regions-manifest-v1-2-bed.zip
fi
