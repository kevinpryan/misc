#!/usr/bin/env nextflow

// stage reference files

/*workflow{
Channel
    .fromPath(params.reference_dir)
    .ifEmpty { error "No reference files found: $params.reference_dir" }
    //.view()
Channel
    .fromFilePairs(params.reads_1)
    .ifEmpty { error "No input files found: $params.reads_1" }
    .set { input_ch }

//input_ch.view()

Channel
    .fromFilePairs(params.reads_2)
    .ifEmpty { error "No input files found: $params.reads_2" }
    .set { input_ch2 }
//input_ch2.view()
// merge the two channels

//input_ch.join(input_ch2, remainder: true)
       //.set { input_ch_merged }
//input_ch.join(input_ch2, by: [0,1]).view()
//input_ch2.view()
input_ch.merge(input_ch2).view()
// .view { it[0] } to view sample names
}
*/

workflow {
    Channel.fromPath("samplesheet_wes_combined.csv")
    | splitCsv( header:true )
    | map { row ->
        meta = id:row.sample
        [meta, [
            file(row.fastq_1, checkIfExists: true),
            file(row.fastq_2, checkIfExists: true)]]
    }
    | view
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
    path reference
    tuple val 

}
*/

