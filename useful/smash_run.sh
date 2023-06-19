#!/bin/bash

#python SMaSH.py -bam /mnt/smash_bams/3532.markdup.sorted.bam,/mnt/smash_bams/3533.markdup.sorted.bam,/mnt/smash_bams/4027.markdup.sorted.bam,/mnt/smash_bams/4028.markdup.sorted.bam -i snps_GRCh38.vcf -o pval_out_test.txt
cd /mnt/SMaSH/
bams=$(ls -m /mnt/data/nf-core-rnaseq-results-subpopulation/star_salmon/*.bam)
python SMaSH.py -bam $bams -i snps_GRCh38.vcf -o pval_out_caf_subpopulation.txt
