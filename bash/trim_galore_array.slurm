#!/bin/bash
#SBATCH --job-name="trim"
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=error/%A_%a.error
#SBATCH -p normal # either normal or highmem - if you need >60GB ram, for for highmem
#SBATCH --mail-user=your_email_address
#SBATCH --mail-type=END
#SBATCH --time=3-11:00:00
#SBATCH --array=1-24 # initiallise 24 jobs - running trim_galore 24 times
#SBATCH --mem-per-cpu=60gb # memory per cpu core


start_time=$SECONDS

cd /data/kryan/sw/NeoFuse

# activate conda environment in which trim_galore is installed
source activate trim

INPUT_DIR=/data/bdigby/Projects/D_O_Connor/raw_reads/
RUN=${SLURM_ARRAY_TASK_ID:-1}


echo "Run: ${RUN}"
read1=$(ls ${INPUT_DIR}*_R1_001.fastq.gz | tail -n 12 | sed -n ${SLURM_ARRAY_TASK_ID}p)
read2=$(ls ${INPUT_DIR}*_R2_001.fastq.gz | tail -n 12 | sed -n ${SLURM_ARRAY_TASK_ID}p)
echo "Read 1: $read1"
echo "Read 2: $read2"

file="$(basename -- $read1)"
echo "Filename: $file"
sample=${file:0:4}
echo "Sample: $sample"
outdir="/data/kryan/rna_seq_bc/trimmed_reads_2nd_half/$sample"
echo "outdir: $outdir"

# run trim galore, can also pass --fastqc to run fastqc once trimming is complete
# comment out this command first to check if your filenames, samples and outdir are as expected using the echo statements above
trim_galore \
	--paired \
	--gzip \
	$read1 \
	$read2 \
	-o $outdir

elapsed=$(( SECONDS - start_time ))
eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"
