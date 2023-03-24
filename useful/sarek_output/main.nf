process MOSDEPTH_OUTPUT_ADD_FILENAME {
    publishDir "$params.outdir/modepth_coverage", mode: 'copy'
    input:
    path(mosdepth_files) // ch_hlatyping

    output:
    path("*.tsv"), emit: ch_hla_samplename_added

    script:
    prefix = mosdepth_files.simpleName 
    """
    awk -v OFS='\t' '{print FILENAME, \$0}' ${mosdepth_files} | tail -n1 >> ${prefix}_with_samplename.tsv
    """
}

ch_mosdepth = Channel.fromPath(mosdepth_files, checkIfExists: true)

MOSDEPTH_OUTPUT_ADD_FILENAME(ch_mosdepth)