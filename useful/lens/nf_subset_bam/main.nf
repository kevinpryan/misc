#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process SORT {
    publishDir "$params.outdir/sorted", mode: 'copy'
    
    input:
    path(bam)

    output:
    path('*.sorted.bam'), emit: ch_bam_sorted

    script:
    prefix = bam.simpleName 
    """
    samtools sort ${bam} -o ${prefix}.sorted.bam
    """
}
process INDEX {
    //label 'samtools_container'

    publishDir "$params.outdir/idx", mode: 'copy'

    input:
    path(bam)

    output:
    //path("*.bam"), emit: ch_bam_indexed
    path("*.bam.bai"), emit: ch_idx

    script:
    """
    samtools index ${bam}
    """
}

process SUBSET_WITH_BED {
    
    //label 'samtools_container'

    publishDir "$params.outdir/subset", mode: 'copy'
    
    input:
    path(bam)
    path(idx)
    path(bed)

    output:
    path("*.subset.bam"), emit: ch_bam_subset
    path("*.subset.bam.bai"), emit: ch_bam_bai_subset

    script:
    prefix = bam.simpleName 
    """
    samtools view -L ${bed} -bX ${bam} ${idx} > ${prefix}.subset.bam
    samtools index ${prefix}.subset.bam
    """
}

// uncomment if you are using the fromPath version of the pipeline
ch_bam = Channel.fromPath(params.bam, checkIfExists: true)
ch_bed = Channel.fromPath(params.bed, checkIfExists: true)

workflow{
    
    SORT(
        ch_bam
    )

    INDEX(
        SORT.out.ch_bam_sorted
    )

    SUBSET_WITH_BED(
        SORT.out.ch_bam_sorted,
        INDEX.out.ch_idx,
        ch_bed
    )   
}
