include { bam2fastq } from '../../../modules/local/bam2fastq'
include { RUN_OPTITYPE } from '../../../modules/local/run_optitype'

workflow optitype{
    /*
    need to realign without alt contigs and remove "chr" from bam header
    */
    take: 
    subsetbam
    reference
    reference_basename

    main:
    bam2fastq(
        subsetbam
    )
    RUN_OPTITYPE(
        bam2fastq.out.subsetfastq,
        dna_rna
    )
    emit:
    RUN_OPTITYPE.out.optitype_call
}
