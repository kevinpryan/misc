#!/bin/bash
#SBATCH --time=3-16:00:00
#SBATCH --job-name="nfcore_rnaseq"
#SBATCH --output=nfcore_rnaseq_inhouse.out
#SBATCH --mail-user=your_email_address
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=20G   # memory per cpu-core
#SBATCH -N 2 # number of nodes

# script to run an nf-core pipeline from a compute node using a specific version of nextflow, using a version of the pipeline
# that I forked from the nf-core github, using a pre-built star index

start_time=$SECONDS

cd /data2/kryan/

echo "Running nf-core/rnaseq on inhouse samples using data2 to hold outputs"

NXF_VER=21.10.6 nextflow run kevinpryan/rnaseq -profile lugh,singularity --input /data2/kryan/inhouse_rnaseq/samplesheet_inhouse_data.csv --fasta /data/kryan/reference/GRCh38.primary_assembly.genome.fa --gtf /data/kryan/reference/gencode.v31.annotation.gtf --max_cpus 12 --outdir /data2/kryan/inhouse_rnaseq/inhouse_data_nfcore_results --deseq2_vst --star_index /data/kryan/rna_seq_bc/caf_subtypes/index/star

elapsed=$(( SECONDS - start_time ))
eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"
