#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process BCFTOOLS_FILTER_MUTECT2_POSITIONS {
    publishDir "${params.outdir}/positions", mode: 'copy'
    input:
    path(vcf)

    output:
    path("*.positions.txt"), emit: ch_positions

    script:
    prefix = vcf.simpleName 
    """
    bcftools query -f '%CHROM\t%POS\n' ${vcf} > ${prefix}.positions.txt
    """
}

process SAMTOOLS_VIEW_FILTER_CRAM_POSITIONS {
    publishDir "${params.outdir}/cram_filtered", mode: 'copy'
    input:
    path(cram)
    path(crai)
    path(positions)
    path(reference)
    path(fai)

    output:
    path("*.bam"), emit: ch_cram_filtered

    script:
    prefix = cram.simpleName 
    """
    samtools view -L ${positions} -b ${cram} --reference ${reference} -t ${fai} > ${prefix}.filtered.mutect2.vars.bam
    """
}

//    samtools view -L ${positions} -b ${cram} -X ${crai} --reference ${reference} > ${prefix}.filtered.mutect2.vars.bam

process SAMTOOLS_SORT {
    publishDir "${params.outdir}/samtools_sorted_index", mode: 'copy'
    input:
    path(bam_filtered)

    output:
    path("*sorted.bam"), emit: ch_bam_sorted

    script:
    prefix = bam_filtered.simpleName
    """
    samtools sort ${bam_filtered} > ${prefix}.sorted.bam
    """ 
}

process SAMTOOLS_INDEX {
    publishDir "${params.outdir}/samtools_sorted_index", mode: 'copy'
    input:
    path(bam_sorted)

    output:
    path("*.bam.bai"), emit: ch_bam_indexed_sorted

    script:
    prefix = bam_sorted.simpleName
    """
    samtools index ${bam_sorted}
    """ 
}


//ch_vcf = Channel.fromPath(params.vcf, checkIfExists: true)
//ch_cram = Channel.fromPath([params.cram1,params.cram2], checkIfExists: true)
ch_vcf = Channel.fromPath(params.vcf, checkIfExists: true)
ch_cram = Channel.fromPath(params.cram1, checkIfExists: true)
ch_cram_crai = Channel.fromPath(params.crai1, checkIfExists: true)
ch_reference = Channel.fromPath(params.reference, checkIfExists: true)
ch_fai = Channel.fromPath(params.reference_idx, checkIfExists: true)

workflow {

    BCFTOOLS_FILTER_MUTECT2_POSITIONS(ch_vcf)
    
    SAMTOOLS_VIEW_FILTER_CRAM_POSITIONS(
        ch_cram,
        ch_cram_crai,
        BCFTOOLS_FILTER_MUTECT2_POSITIONS.out.ch_positions,
        ch_reference,
        ch_fai
    )

    SAMTOOLS_SORT(
        SAMTOOLS_VIEW_FILTER_CRAM_POSITIONS.out.ch_cram_filtered
    )

    SAMTOOLS_INDEX(
       SAMTOOLS_SORT.out.ch_bam_sorted
    )
}


