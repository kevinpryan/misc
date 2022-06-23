#!/bin/bash
#SBATCH --job-name="strandedness"
#SBATCH --output=/data/kryan/rna_seq_bc/batch/logs/%A_%a.out
#SBATCH --error=/data/kryan/rna_seq_bc/batch/error/%A_%a.error
#SBATCH -p normal
#SBATCH --mail-user=k.ryan45@nuigalway.ie
#SBATCH --mail-type=END
#SBATCH --time=3-16:00:00
#SBATCH --mem-per-cpu=60gb # memory per cpu core
#SBATCH --array=1-34 # create a 28-task job array

# script to run strandedness
source activate strandedness

cd /data/kryan/rna_seq_bc/caf_subtypes/EGAD00001004810/strandedness_test

echo "Argument provided: $1"

sample=$(cut -f1 $1 | tail -n +2 | sed -n ${SLURM_ARRAY_TASK_ID}p)
read1=$(cut -f2 $1 | tail -n +2 | sed -n ${SLURM_ARRAY_TASK_ID}p)
read2=$(cut -f3 $1 | tail -n +2 | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo "Sample: $sample"
echo "Read 1: $read1"
echo "Read 2: $read2"

check_strandedness -g /data/kryan/reference/gencode.v31.annotation.gtf -fa /data/kryan/reference/gencode.v40.transcripts.fa -r1 $read1 -r2 $read2 -k /data/kryan/rna_seq_bc/caf_subtypes/EGAD00001005744/kallisto_index > ${sample}_strandedness.txt


