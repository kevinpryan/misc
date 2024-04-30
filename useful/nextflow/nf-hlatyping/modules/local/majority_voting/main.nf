process MAJORITY_VOTE{
    publishDir "$params.outdir/results"

    input:
    tuple val(meta), path(outputs)
    path benchmark

    output:
    tuple val(meta), path("*_all_calls_mhci.tsv"), emit: all_calls
    tuple val(meta), path("*_majority_vote_mhci.tsv"), emit: majority_vote

    script:
    """
    Rscript parse_outputs_majority_vote.R --samplename ${meta.sample} --optitype ${meta.sample} --polysolver ${meta.sample} --hlala ${meta.sample}/${meta.sample}/hla/ --kourami ${meta.sample} --benchmark ${benchmark}
    """
}