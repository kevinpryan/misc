#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GET_NEW_NAMES {
    publishDir "$params.outdir/new_names", mode: 'copy'

    input:
    path(filedir)
    output:
    path("samplesheet.csv"), emit: ch_samplesheet


    script:
    """
    # can swap this script out for another script depending on your renaming needs
    Rscript '${baseDir}/bin/rename_files_general.R' -f $filedir 
    """

}

process GENERATE_MANIFEST {
    publishDir "$params.outdir/manifest", mode: 'copy'

    input:
    path(samplesheet)
    val(file_out)
    val(dataset_name)

    output:
    path("*.csv"), emit: ch_manifest

    script:
    """
    python '${baseDir}/bin/generate_manifest.py' $samplesheet $file_out -d $dataset_name
    """
}
process RENAME_FILES{
    publishDir "$params.outdir/copied_files", mode: 'copy'
    input:
    //path(inpath)
    path(new_filenames) // output of rename_files_general.R - GET_NEW_NAMES
    val(dna_seq_type) // params.dna_seq_method
    val(rna_seq_type) // params.rna_seq_method

    output:
    path("*.fastq.gz"), optional: true

    shell:
    """
    #!/bin/bash
    bash '${baseDir}/bin/rename.sh' -f $new_filenames -d $dna_seq_type -r $rna_seq_type
     """
}


// uncomment if you are using the fromPath version of the pipeline
ch_inpath = Channel.fromPath(params.input_dir, checkIfExists: true)
ch_infile = Channel.of(params.samplesheet)
ch_manifest_name = Channel.of(params.manifest_name)
ch_dataset_name = Channel.of(params.dataset_name)
ch_input_dir_dna = Channel.fromPath(params.input_dir_dna, checkIfExists: true)
ch_input_dir_rna = Channel.fromPath(params.input_dir_rna, checkIfExists: true)

workflow{

    GET_NEW_NAMES(
        ch_infile
    )

    GENERATE_MANIFEST(
        GET_NEW_NAMES.out.ch_samplesheet,
        ch_manifest_name,
        ch_dataset_name
    )

     RENAME_FILES(
        GET_NEW_NAMES.out.ch_samplesheet,
        ch_input_dir_dna,
        ch_input_dir_rna
    )
    
    
}
