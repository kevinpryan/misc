#!/usr/bin/env nextflow

include {alt_align} from "../subworkflows/local/alt_align"
include { prepPolysolver } from "../subworkflows/local/prepPolysolver"
include { optitype } from "../subworkflows/local/optitype"
include { polysolver } from "../subworkflows/local/polysolver"
include { hlala } from "../subworkflows/local/hlala"
include { kourami } from "../subworkflows/local/kourami"
include { FASTP } from "../modules/nf-core/fastp"

workflow HLATYPING {
    // TODO: add samplesheet check, seq_type should be in dna,rna 
    take:
    ch_fastq
    reference_basename
    reference_dir
    hlatypes
    chr
    hla_la_graph
    kourami_ref
    kourami_database
    trimmer
    save_trimmed_fail
    save_merged
    adapter_fasta
    
    main:
    reference_basename = Channel.value(reference_basename)
    ch_ref = file(reference_dir, checkIfExists: true)
    ch_hlatypes = file(hlatypes, checkIfExists: true)
    chromosome = Channel.value(chr)
    ch_graph = file(hla_la_graph, checkIfExists: true)
    ch_ref_kourami = file(kourami_ref, checkIfExists: true)
    ch_db_kourami = file(kourami_database, checkIfExists: true)
    ch_fastq.view()
    if (trimmer == 'fastp') {
    //ch_adapter_fasta = Channel.empty()
    FASTP (
    ch_fastq,
    [],
    save_trimmed_fail,
    save_merged
    )
    ch_fastq_align = FASTP.out.reads
    } else {
    ch_fastq_align = ch_fastq
    }
    alt_align(
    ch_fastq_align,
    ch_ref,
    ch_hlatypes,
    reference_basename,
    chromosome
    )
    optitype(
        alt_align.out//,
        //dna_rna
    )
    polysolver(
        alt_align.out,
        ch_ref,
        reference_basename
    )
    hlala(
        alt_align.out,
        ch_graph
    ) 
    kourami(
        alt_align.out,
        ch_db_kourami,
        ch_ref_kourami
    )
    // rough idea of end of pipeline
    // this will be in a subworkflow
    // see /home/kevin/Documents/PhD/nextflow_test/test2/main.nf for working example
    // not sure whether to use arcashla or not
    // docker image: r-basic:dev
    RUN_OPTITYPE.out.mix(RUN_KOURAMI.out, RUN_POLYSOLVER.out, RUN_HLALA.out)
           .groupTuple(by: 0, size: 3)
           .set(HLATYPING_RESULTS)
           //.view()
    COMBINE_HLATYPING(HLATYPING_RESULTS)
    ARCASHLA_CONVERT(COMBINE_HLATYPING.out.table)
    MAJORITY_VOTE(ARCASHLA_CONVERT.out.converted)
}
