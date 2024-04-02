#!/usr/bin/env nextflow

include {alt_align} from "./workflows/local/alt_align"
include { prepPolysolver } from "./workflows/local/prepPolysolver"

/*
process bwa_mem_align_alt{
    publishDir "$params.outdir/align"

    input:
    path reference
    tuple val(meta), path(reads)
    val reference_basename

    output:
    tuple val(meta), path("*.bwamem.sam"), emit: samfile
    path("*.bwamem.sam.header")
    path("*.bwamem.sam.flagstat")
    script:
    """
    bwa mem -t ${task.cpus} ${reference_basename}.fa ${reads} > "${meta.sample}.bwamem.sam"
    samtools view -H ${meta.sample}.bwamem.sam > ${meta.sample}.bwamem.sam.header
    samtools flagstat ${meta.sample}.bwamem.sam > ${meta.sample}.bwamem.sam.flagstat
    """
}
*/

///*
//process bwa_mem_postalt{
//    publishDir "$params.outdir/postalt"
//
//    input:
//    path reference
//    tuple val(meta), path(samfile)
//    val reference_basename
//    output:
//    tuple val(meta), path("*_postalt.bam"), emit: bamfile_postalt
//    path("*_postalt.sam.header")
//    path("*_postalt.sam.flagstat")
//    path("*_postalt.bam.flagstat")
//
//    script:
//    """
//    k8-linux /usr/local/bin/bwa-0.7.15/bwa-postalt.js \
//    ${reference_basename}.fa.alt \
//    ${meta.sample}.bwamem.sam > ${meta.sample}_postalt.sam
//    samtools view -H ${meta.sample}_postalt.sam > ${meta.sample}_postalt.sam.header
//    samtools flagstat ${meta.sample}_postalt.sam > ${meta.sample}_postalt.sam.flagstat
//    samtools view -bh -o ${meta.sample}_postalt.bam ${meta.sample}_postalt.sam
//    samtools flagstat ${meta.sample}_postalt.bam > ${meta.sample}_postalt.bam.flagstat
//    """
//}
//*/

/*
process bam2fastq{
    publishDir "$params.outdir/bam2fastq"

    input:
    tuple val(meta), path(subsetbam), path(subsetbam_bai)
    output:
    tuple val(meta), path("*.fq"), emit: subsetfastq

    script:
    """
    sambamba view \
    -f "bam" -h -p -l 0 -t ${task.cpus} \
    ${subsetbam} | sambamba sort -p -n -t ${task.cpus} -o - /dev/stdin |
    samtools fastq /dev/stdin \
    -1 "${meta.sample}_subset.1.fq" \
    -2 "${meta.sample}_subset.2.fq" \
    -0 /dev/null -s /dev/null -n
    """
}
*/

/*
process realignwithoutAlt{
    publishDir "$params.outdir/realignwithoutAlt"
    
    input:
    tuple val(meta), path(reads)
    path reference
    val reference_basename

    output:
    tuple val(meta), path("*_realign.sam"), emit: realignbam
    path("*_realign.sam.header")
    path("*_realign.sam.flagstat")

    script:
    """
    bwa mem -t ${task.cpus} -j ${reference_basename}.fa ${reads} > ${meta.sample}_realign.sam
    samtools view -H ${meta.sample}_realign.sam > ${meta.sample}_realign.sam.header
    samtools flagstat ${meta.sample}_realign.sam > ${meta.sample}_realign.sam.flagstat
    """
}
*/

/*
process reheaderChr{
    publishDir "$params.outdir/reheaderChr"

    input:
    tuple val(meta), path(bamfile), path(bamfile_idx)

    output:
    //tuple val(meta), path("*_reheader.bam"), path("*_reheader.bam.bai"), path("*_reheader.bam.flagstat"), path("*_reheader.bam.header"), emit: reheaderbam
    tuple val(meta), path("*_reheader.bam"), path("*_reheader.bam.bai"), emit: reheaderbam
    path("*_reheader.bam.flagstat")
    path("*_reheader.bam.header")
    script:
    """
    bash reheader_chr.sh ${bamfile} ${meta.sample}_reheader.bam
    samtools flagstat ${meta.sample}_reheader.bam > ${meta.sample}_reheader.bam.flagstat
    samtools view -H ${meta.sample}_reheader.bam > ${meta.sample}_reheader.bam.header
    """
}
*/

process RUN_OPTITYPE{
    publishDir "$params.outdir/OPTITYPE/${meta.sample}"

    input:
    tuple val(meta), path(subsetfastq)
    val dna_rna

    output:
    tuple val(meta), path("*"), emit: optitype

    script:
    """
    OptiTypePipeline.py --input *1.fq *.2.fq --verbose --${dna_rna} --outdir /data4/kryan/misc/useful/nextflow/nf-hlatyping/testdir/3532

    optitype -i ${subsetfastq} --${dna_rna} --outdir ./
    """
}

workflow OPTITYPE{

    take: 
    subsetbam
    reference
    reference_basename
    dna_rna

    main:
    bam2fastq(
        subsetbam
    )
    RUN_OPTITYPE(
        bam2fastq.out.subsetfastq,
        dna_rna
    )

}

// this version of the workflow works 20240329

workflow {
    Channel.fromPath(params.samplesheet, checkIfExists: true)
    | splitCsv( header:true )
    | map { row ->
        meta = row.subMap('sample')
        [meta, [
            file(row.fastq_1, checkIfExists: true),
            file(row.fastq_2, checkIfExists: true)]]
    }
    | set { ch_fastq }
    reference_basename = Channel.value(params.reference_basename)
    ch_ref = file(params.reference_dir, checkIfExists: true)
    ch_hlatypes = file(params.hlatypes, checkIfExists: true)
    chromosome = Channel.value(params.chr)
    dna_rna = Channel.value(params.dna_rna)
    alt_align(
    ch_fastq,
    ch_ref,
    ch_hlatypes,
    reference_basename,
    chromosome
    )
    prepPolysolver(
    alt_align.out,
    ch_ref,
    reference_basename
    )
}


