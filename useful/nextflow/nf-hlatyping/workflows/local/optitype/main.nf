include { bam2fastq } from '../../../modules/local/bam2fastq'
include { RUN_OPTITYPE } from '../../../modules/local/run_optitype'

workflow optitype{
    /*
    convert bam to fastq then run Optitype
    */
    take: 
    subsetbam
    dna_rna

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
