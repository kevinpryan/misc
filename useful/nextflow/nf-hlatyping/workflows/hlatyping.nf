#!/usr/bin/env nextflow

include {alt_align} from "../subworkflows/local/alt_align"
include { prepPolysolver } from "../subworkflows/local/prepPolysolver"
include { optitype } from "../subworkflows/local/optitype"
include { polysolver } from "../subworkflows/local/polysolver"
include { hlala } from "../subworkflows/local/hlala"
include { kourami } from "../subworkflows/local/kourami"
include { FASTP } from "../modules/nf-core/fastp"
include { MAJORITY_VOTE } from "../modules/nf-core/fastp"

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
    benchmarking
    
    main:
    reference_basename = Channel.value(reference_basename)
    ch_ref = file(reference_dir, checkIfExists: true)
    ch_hlatypes = file(hlatypes, checkIfExists: true)
    chromosome = Channel.value(chr)
    ch_graph = file(hla_la_graph, checkIfExists: true)
    ch_ref_kourami = file(kourami_ref, checkIfExists: true)
    ch_db_kourami = file(kourami_database, checkIfExists: true)
    ch_benchmark = file("$projectDir/assets/benchmarking_results_claeys.csv", checkIfExists: true)
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

    // untested on real data from here on in
    RUN_OPTITYPE.out.optitype.mix(RUN_KOURAMI.out.kourami, RUN_POLYSOLVER.out.polysolver)
           .groupTuple(by: 0, size: 3)
           .set{ ch_hlatyping_outputs }
    ch_hlatyping_outputs
                    .map{meta, results ->
                        [ meta, results.collect { it.getParent() } ]
                    }
                    .set{ ch_hlatyping_outputs_grouped }
    MAJORITY_VOTE(
        ch_hlatyping_outputs_grouped,
        ch_benchmark
    )
    MAJORITY_VOTE.out.majority_vote.collectFile(storeDir: "${params.outdir}/combined_results", name: 'nf_core_hlatyping_results_majority_vote.tsv', newLine: true, keepHeader: 1, sort: { it[0] }) { it[1] }
    MAJORITY_VOTE.out.all_calls.collectFile(storeDir: "${params.outdir}/combined_results", name: 'nf_core_hlatyping_results_all_calls.tsv', newLine: true, keepHeader: 1, sort: { it[0] }) { it[1] }   
    
}
