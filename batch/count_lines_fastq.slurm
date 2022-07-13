#!/bin/bash
#SBATCH --time=3-16:00:00
#SBATCH --job-name="count_fastq"
#SBATCH --output=fastq_countlines.out
#SBATCH --mail-user=your_email_address
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=60G   # memory per cpu-core

# date written = 14/02/2022
# script to count number of lines in gzipped fastq files in a particular directory

start_time=$SECONDS

for x in /data/bdigby/Projects/D_O_Connor/raw_reads/*.fastq.gz; do zcat "$x" | wc -l && echo $x; done

eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"
