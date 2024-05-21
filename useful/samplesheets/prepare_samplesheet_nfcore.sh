#!/bin/bash
conda activate nf-core
set -e
cd /data4/caf_rna_bc/lung_cancer_novogene/rna_seq
python fastq_dir_to_samplesheet.py /data4/caf_rna_bc/lung_cancer_novogene/rna_seq/X204SC23110279-Z01-F001_01 samplesheet_rna_lung1_nfcore.csv -re -st reverse --read1_extension "1.fq.gz" --read2_extension "2.fq.gz" -sn -si 2
python fastq_dir_to_samplesheet.py /data4/caf_rna_bc/lung_cancer_novogene/rna_seq/X204SC23110279-Z01-F001_02 samplesheet_rna_lung2_nfcore.csv -re -st reverse --read1_extension "1.fq.gz" --read2_extension "2.fq.gz" -sn -si 2
tail -n +2 samplesheet_rna_lung2_nfcore.csv > samplesheet_rna_lung2_nfcore_remove_header.csv
cat samplesheet_rna_lung1_nfcore.csv samplesheet_rna_lung2_nfcore_remove_header.csv > samplesheet_rna_combined_nfcore.csv
#cat samplesheet_rna_combined.csv | cut -f 1-3 -d "," > samplesheet_rna_combined_remove_column.csv
rm samplesheet_rna_lung2_nfcore_remove_header.csv samplesheet_rna_lung2_nfcore.csv samplesheet_rna_lung1_nfcore.csv
