#!/bin/bash
#SBATCH --time=16:00:00
#SBATCH --job-name="fisher_cov_w"
#SBATCH --output=fisher_enrichment_whole_blood_nocov.out
cd /data/kryan/project/gtex/scripts/

source activate renv4_1

Rscript --vanilla enrichment_analysis.R "default" "whole_blood_norm" "pass"
