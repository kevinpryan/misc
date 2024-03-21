#!/bin/bash

start_time=$SECONDS
# Check if the correct number of arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 input_directory output_directory reference_file"
    exit 1
fi

# Assign input arguments to variables
indir="$1"
outdir="$2"
reference="$3"

set -o pipefail
set -e
module load Anaconda3/2024.02-1
conda activate bwa
mkdir -p "$outdir"
cd "$outdir"
echo "current working directory..."
pwd
# For 'base', to match files ending with either fastq.gz or altalt.read1.fq
base=$(find "${indir}" -type f \( -name '*.fastq.gz' -o -name 'altalt.read1.fq' \) -exec basename {} \; | cut -f 1 -d "_" | head -n 1)

# For 'fq1', to match files containing '*R1*fastq.gz' or 'altalt.read1.fq'
fq1=$(find "${indir}" -type f \( -name '*R1*fastq.gz' -o -name 'altalt.read1.fq' \) -exec basename {} \; | head -n 1)

# For 'fq2', to match files containing '*R2*fastq.gz' or 'altalt.read2.fq'
fq2=$(find "${indir}" -type f \( -name '*R2*fastq.gz' -o -name 'altalt.read2.fq' \) -exec basename {} \; | head -n 1)

#base=$(echo "altalt.read1.fq" | cut -f 1 -d ".")
echo "samplename..."
echo $base

echo "fq1..."
echo $fq1

echo "fq2..."
echo $fq2

outdir_sample="${outdir}/${base}"
basename_out="${base}_postalt"

mkdir -p "${outdir_sample}"
export PATH=$PATH:/data/kryan/bwakit-0.7.15/bwa-0.7.15/
mkdir -p "${outdir_sample}/bwamem"

echo "aligning..."
bwa mem -t 15 "${reference}" "${indir}/${fq1}" "${indir}/${fq2}" > "${outdir_sample}/bwamem/${base}.bwamem.sam"
samtools view -H ${outdir_sample}/bwamem/${base}.bwamem.sam > ${outdir_sample}/bwamem/${base}.bwamem.sam.header
samtools flagstat ${outdir_sample}/bwamem/${base}.bwamem.sam > ${outdir_sample}/bwamem/${base}.bwamem.sam.flagstat
export PATH=$PATH:/data/kryan/bwakit-0.7.15/k8
export PATH=$PATH:/data/kryan/bwakit-0.7.15/bwa-0.7.15/bwakit

mkdir -p "${outdir_sample}/postalt"
echo "postalt..."
k8-linux /data/kryan/bwakit-0.7.15/bwa-0.7.15/bwakit/bwa-postalt.js \
    ${reference}.alt \
    ${outdir_sample}/bwamem/${base}.bwamem.sam > ${outdir_sample}/postalt/${basename_out}.sam
samtools view -H ${outdir_sample}/postalt/${basename_out}.sam > ${outdir_sample}/postalt/${basename_out}.sam.header
samtools flagstat ${outdir_sample}/postalt/${basename_out}.sam > ${outdir_sample}/postalt/${basename_out}.sam.flagstat

rm "${outdir_sample}/bwamem/${base}.bwamem.sam"

mkdir -p "${outdir_sample}/samtools_view"
samtools view -bh -o "${outdir_sample}/samtools_view/${basename_out}.bam" "${outdir_sample}/postalt/${basename_out}.sam"
samtools flagstat "${outdir_sample}/samtools_view/${basename_out}.bam" > "${outdir_sample}/samtools_view/${basename_out}.bam.flagstat"
echo "Post-alignment finished"

rm "${outdir_sample}/postalt/${basename_out}.sam"

sambamba sort -o "${outdir_sample}/samtools_view/${basename_out}_sorted.bam" "${outdir_sample}/samtools_view/${basename_out}.bam"
samtools flagstat "${outdir_sample}/samtools_view/${basename_out}_sorted.bam" > "${outdir_sample}/samtools_view/${basename_out}_sorted.bam.flagstat"
samtools index "${outdir_sample}/samtools_view/${basename_out}_sorted.bam"

bammarkduplicates I="${outdir_sample}/samtools_view/${basename_out}_sorted.bam" O="${outdir_sample}/samtools_view/${basename_out}_sorted_mdup.bam" index=1 rmdup=0

echo "subsetting BAM"
mkdir -p "${outdir_sample}/subset-alt"
grep 'chr6_' ${reference}.alt | awk -F '\t' '{print $1}' > "${outdir_sample}/subset-alt/bwakit-alt_chr6_contigs.txt"
awk -F '\t' '{print $3}' /data/kryan/bwakit-0.7.15/bwa-0.7.15/bwakit/bwa.kit/resource-human-HLA/HLA-ALT-type.txt > "${outdir_sample}/subset-alt/bwakit-hla_contigs.txt"
echo "contents of subset-alt dir before subsetting..."
ls -l "${outdir_sample}/subset-alt/"

primary_contig='chr6:28509970-33480727'
alt_contigs=$(<"${outdir_sample}/subset-alt/bwakit-alt_chr6_contigs.txt")
hla_contigs=$(<"${outdir_sample}/subset-alt/bwakit-hla_contigs.txt")
no_contig='*'
OUT="${outdir_sample}/subset-alt/${basename_out}_subset_mdup.bam"
OUT_SORTED="${outdir_sample}/subset-alt/${basename_out}_subset_mdup_sorted.bam"
echo "OUT = "
echo $OUT
# subset the bam in which duplicates have been marked to only include contigs from chromosome 6, alt chromosome 6 contigs or unmapped contigs
samtools view -o "${OUT}" -b "${outdir_sample}/samtools_view/${basename_out}_sorted_mdup.bam" $primary_contig $alt_contigs $hla_contigs "$no_contig" 
samtools sort -o "${OUT_SORTED}" "${OUT}"
sambamba index "${OUT_SORTED}"
samtools view -H ${OUT_SORTED} > ${OUT_SORTED}.header
samtools flagstat ${OUT_SORTED} > ${OUT_SORTED}.flagstat

elapsed=$(( SECONDS - start_time ))
eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"
