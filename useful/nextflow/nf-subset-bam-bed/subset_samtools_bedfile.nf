#!/usr/bin/env nextflow

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


process SUBSET {
    publishDir "$params.outdir/chr21", mode: 'copy'
    
    input:
    path(bamfile)
    path(idx)

    output:
    path("*.chr21.bam"), emit: subsample

    script:
    prefix = bamfile.simpleName
    """
    samtools view -b -h -L ${bed} ${bam} > ${prefix}.subset.bam
    """    
}


ch_bams = Channel.fromPath(['/data4/kryan/rnafusion_inhouse_aug2023/outputs_pt7-9/star_for_arriba/3532.Aligned.out.bam',
'/data4/kryan/rnafusion_inhouse_aug2023/outputs_pt4-6/star_for_arriba/4315.Aligned.out.bam',
'/data4/kryan/rnafusion_inhouse_aug2023/outputs_pt7-9/star_for_arriba/4340.Aligned.out.bam',
'/data4/kryan/rnafusion_inhouse_aug2023/outputs_pt7-9/star_for_arriba/4344.Aligned.out.bam'],
checkIfExists: true)

ch_bams_map = ch_bams.map { it -> [it.simpleName.toString(), it]}

ch_beds = Channel.fromPath(['/data4/kryan/misc/useful/R/3532/3532.bed',
'/data4/kryan/misc/useful/R/4315/4315.bed',
'/data4/kryan/misc/useful/R/4340/4340.bed',
'/data4/kryan/misc/useful/R/4344/4344.bed'],
checkIfExists: true)

ch_beds_map = ch_beds.map { it -> [it.simpleName.toString(), it ] }
workflow {
/*
    INDEX(ch_bams)
    SUBSET(
    ch_bams,
    INDEX.out.ch_idx
    )
*/
}

