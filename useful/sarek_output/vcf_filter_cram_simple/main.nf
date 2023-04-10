#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process BCFTOOLS_FILTER_MUTECT2_POSITIONS {
    publishDir "$params.outdir/positions", mode: 'copy'
    input:
    path(vcf)

    output:
    path(positions), emit: ch_positions

    script:
    prefix = vcf.simpleName 
    """
    bcftools query -f '%CHROM\t%POS\n' ${vcf} > ${prefix}.positions.txt
    """
}

process SAMTOOLS_VIEW_FILTER_CRAM_POSITIONS {
    publishDir "$params.outdir/cram_filtered", mode: 'copy'
    input:
    path(cram)
    path(positions)

    output:
    path(*.bam), emit: ch_cram_filtered

    script:
    prefix = cram.simpleName 
    """
    samtools view -L ${positions} -b ${cram} > ${prefix}.filtered.mutect2.vars.bam
    """
}

process SAMTOOLS_SORT {
    publishDir "$params.outdir/samtools_sort", mode: 'copy'
    input:
    path(bam_filtered)

    output:
    path(*.bam), emit: ch_bam_sorted

    script:
    prefix = bam_filtered.simpleName
    """
    samtools sort ${bam_filtered}
    """ 
}

process SAMTOOLS_INDEX {
    publishDir "$params.outdir/samtools_index", mode: 'copy'
    input:
    path(bam_sorted)

    output:
    path(*.bam.*), emit: ch_bam_indexed_sorted

    script:
    prefix = bam_sorted.simpleName
    """
    samtools index ${cram_filtered}
    """ 
}

ch_vcf = Channel.fromPath(params.vcf, checkIfExists: true)
ch_cram = Channel.fromPath(params.cram, checkIfExists: true)

workflow {

    BCFTOOLS_FILTER_MUTECT2_POSITIONS(ch_vcf)
    SAMTOOLS_VIEW_FILTER_CRAM_POSITIONS(
        ch_cram,
        BCFTOOLS_FILTER_MUTECT2_POSTIONS.out.ch_positions
    )
    SAMTOOLS_SORT(
        SAMTOOLS_VIEW_FILTER_CRAM_POSITIONS.out.ch_cram_filtered
    )
    SAMTOOLS_INDEX(
        SAMTOOLS_SORT.out.ch_bam_indexed_sorted
    )
}


