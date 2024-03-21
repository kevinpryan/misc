# get CTA genes used in LENS
wget https://gitlab.com/landscape-of-effective-neoantigens-software/nextflow/modules/tools/lens/-/wikis/uploads/5a9786203497b90c0cc0c0a6a251399b/cta_and_self_antigen.homo_sapiens.gene_list
# get BED file for the enrichment panel
wget https://emea.support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/illumina-prep/exome/hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED 
# count number of genes in CTA list 
cat cta_and_self_antigen.homo_sapiens.gene_list | sort | uniq | wc -l

