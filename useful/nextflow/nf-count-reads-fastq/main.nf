#!/usr/bin/env nextflow
params.reads = '/home/kevin/Documents/PhD/wes_bc/dummy_files/subset_rnaseq_12lines/sub_ar-sub*_{1,2}.fastq.gz.gz'
//params.reads = '/home/kevin/Documents/PhD/wes_bc/dummy_files/dummy_files_all/*_R{1,2}_001.fastq.gz'
params.outdir = "."
process COUNT_READS{
    publishDir "$params.outdir/fastq_count", mode: 'copy'
    input:
    tuple val(sample), path(fastq)
    output:
    path("*.txt"), emit: ch_txt

    script:
    """
    echo "sample: ${sample} read 1" > "${sample}.txt"
    zcat ${fastq[0]} | echo \$((`wc -l`/4)) >> "${sample}.txt"
    echo "read 2" >> "${sample}.txt"
    zcat ${fastq[1]} | echo \$((`wc -l`/4)) >> "${sample}.txt"
    """
}

workflow{
Channel
  .fromFilePairs( params.reads, checkIfExists: true )
  .set { read_pairs_ch }


COUNT_READS(read_pairs_ch)
//read_pairs_ch.view{ it[0] } - view sample id 
//read_pairs_ch.view{ it[1][0]} - view read 1
//read_pairs_ch.view{ it[1][1]} - view read 2
}