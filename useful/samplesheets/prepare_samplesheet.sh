#!/bin/bash
cd /data4/caf_rna_bc/lung_cancer_novogene/data_download
python fastq_dir_to_samplesheet.py /data4/caf_rna_bc/lung_cancer_novogene/data_download/X204SC23110279-Z01-F001_01 samplesheet_rna_lung1.csv -re --read1_extension 1.fq.gz --read2_extension 2.fq.gz -sn -si 2
python fastq_dir_to_samplesheet.py /data4/caf_rna_bc/lung_cancer_novogene/data_download/X204SC23110279-Z01-F001_02 samplesheet_rna_lung2.csv -re --read1_extension 1.fq.gz --read2_extension 2.fq.gz -sn -si 2
tail -n +2 samplesheet_rna_lung2.csv > samplesheet_rna_lung2_remove_header.csv
cat samplesheet_rna_lung1.csv samplesheet_rna_lung2_remove_header.csv > samplesheet_rna_combined.csv
cat samplesheet_rna_combined.csv | cut -f 1-3 -d "," > samplesheet_rna_combined_remove_column.csv
rm samplesheet_rna_combined.csv samplesheet_rna_lung1.csv samplesheet_rna_lung2.csv samplesheet_rna_lung2_remove_header.csv
