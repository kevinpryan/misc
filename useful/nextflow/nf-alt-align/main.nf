#!/usr/bin/env nextflow
params.reads_1 = '/data2/caf_wes_bc/raw_data/*{R1,R2}*.fastq.gz'
params.reads_2 = '/data2/caf_wes_bc/raw_data_final_samples/*{R1,R2}*.fastq.gz'

Channel
    .fromFilePairs(params.reads_1, flat: true)
    .ifEmpty { error "No input files found: $params.reads_1" }
    .set { input_ch }

//input_ch.view()

Channel
    .fromFilePairs(params.reads_2, flat: true)
    .ifEmpty { error "No input files found: $params.reads_2" }
    .set { input_ch2 }
//input_ch2.view()
// merge the two channels
input_ch.join(input_ch2, remainder: true).view { it[0] }
    //.set { input_ch_merged }

//input_ch_merged.view() 

process INDEX {
    publishDir "$params.outdir/idx", mode: 'copy'

    input:
    path(bamfile)

    output:
    path("*.bam.bai"), emit: ch_idx

    script:
    """
    samtools index ${bamfile}
    """
}

process bwa_mem_align_alt{
    publishDir "$params.outdir/align"

    input:
    path()

}