#!/usr/bin/env nextflow

// stage reference files
workflow{
Channel
    .fromPath(params.reference_dir)
    .ifEmpty { error "No reference files found: $params.reference_dir" }
    .view()
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
input_ch.join(input_ch2, remainder: true)
    .set { input_ch_merged }
// .view { it[0] } to view sample names
}
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

//bwa mem -t 15 "${reference}" "${indir}/${fq1}" "${indir}/${fq2}" > "${outdir_sample}/bwamem/${base}.bwamem.sam"
///samtools view -H ${outdir_sample}/bwamem/${base}.bwamem.sam > ${outdir_sample}/bwamem/${base}.bwamem.sam.header
//samtools flagstat ${outdir_sample}/bwamem/${base}.bwamem.sam > ${outdir_sample}/bwamem/${base}.bwamem.sam.flagstat

/*
process bwa_mem_align_alt{
    publishDir "$params.outdir/align"

    input:
    path()

}
*/
