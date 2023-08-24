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
    samtools view -bX ${bamfile} ${idx} chr21 > ${prefix}.chr21.bam
    """    
}

ch_bams = Channel.fromPath(['/data2/kryan/lens/raft_install/projects/retry-my-lens-proj/work/46/2fb3dc54adecce0302ebf503165cfd/caf_neo_bc-Pt-06-ar-4315.sorted.bam',
'/data2/kryan/lens/raft_install/projects/retry-my-lens-proj/work/84/d46ce665ae576533d5775999eadbf5/caf_neo_bc-Pt-06-nd-4316.sorted.bam',
'/data2/kryan/lens/raft_install/projects/retry-my-lens-proj/work/88/4e72355c658ae8dac7dffdef492e40/caf_neo_bc-Pt-06-nr-4316.Aligned.sortedByCoord.out.bam'],
checkIfExists: true)

workflow {
    INDEX(ch_bams)
    SUBSET(
    ch_bams,
    INDEX.out.ch_idx
    )
}

