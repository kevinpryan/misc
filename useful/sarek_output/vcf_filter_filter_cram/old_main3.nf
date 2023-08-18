#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Read in the list of VCF files

//vcf_files = Channel.fromPath(params.vcf_path)
//    .map { it -> [ name: it.simpleName, fullpath: it, tan: it.simpleName.toString().split('\\.')[0].split('\\_')[0], caf: it.simpleName.toString().split('\\.')[0].split('\\_')[2] ] }
vcf_files_caf = Channel.fromPath(params.vcf_path)
    //.map { it -> [ it.simpleName.toString().split('\\.')[0].split('\\_')[0], fullpath: it ] }
    .map { it -> [ it.simpleName.toString().split('\\.')[0].split('\\_')[0], it, it.simpleName ] }

vcf_files_tan = Channel.fromPath(params.vcf_path)
    //.map { it -> [ it.simpleName.toString().split('\\.')[0].split('\\_')[2], fullpath: it ] }
    .map { it -> [ it.simpleName.toString().split('\\.')[0].split('\\_')[2], it, it.simpleName ] }
    //.view()
// Read in the list of CRAM files
cram_files = Channel.fromPath(params.cram_path)
     .map{ it -> [ it.simpleName.toString(), it]}
     //.view()
crai_files = Channel.fromPath(params.crai_path)
//select     .map{ it -> [ it.simpleName.toString(), it]}

cram_crai = cram_files
            .join(crai_files)
            //.view()
// Combine the VCF and CRAM files based on the sample IDs
vcf_cram_pairs_caf = cram_files
                     .join(vcf_files_caf)
                     //.view()

vcf_cram_pairs_tan = cram_files
                     .join(vcf_files_tan)
                     //.view()

vcf_cram_combined_full = vcf_cram_pairs_caf
                         .join(vcf_cram_pairs_tan, by:3)
                         .unique()
                         //.view()
//norm_vcfs_input = vcf_cram_pairs_caf
 //                 .join(vcf_cram_pairs_tan, by:3)
//                  .collect { it[3], it[2] }
//                  .view()

norm_vcfs_input = vcf_cram_pairs_caf
                  .map { it -> [ it[3], it[2]]}
                  .groupTuple()
                  //.collect()
                  //.view()
ch_fasta = Channel.fromPath(params.fasta, checkIfExists: true)
                  .collect()


process NORMALISE_VCFS {
    label 'bcftools_container'
    publishDir "$params.outdir/norm_vcf", mode: 'copy'
    input:
    tuple val(patient), path(vcf)
    path(fasta)
    output:
    tuple val(patient), path("*_normalised.vcf.gz"), emit: vcf_norm
    script:
    """
    bcftools norm --fasta-ref $fasta --output-type z --output ${vcf.simpleName}_normalised.vcf.gz $vcf 
    """    
}


process VCF2BED {
    label 'bedops_container'
    publishDir "$params.outdir/vcf2bed", mode: 'copy'
    input:
    tuple val(patient), path(vcf_norm)
    output:
    tuple val(patient), path("*.bed"), emit: bed_from_vcf
    script:
    """
    vcf2bed < <( gunzip -c $vcf_norm ) | cut -f 1-3 > ${patient}.bed 
    """
}


workflow {
    //NORMALISE_VCFS(norm_vcfs_input,
    //               ch_fasta)
    // vcf_norm is [1], cram is [3] in this channel
    //vcf_cram_combined_full_norm = NORMALISE_VCFS.out.vcf_norm.join(vcf_cram_combined_full).map {it -> [it[0], it[1], it[3]] }
//											  .groupTuple()
//VCF2BED(NORMALISE_VCFS.out.vcf_norm)
VCF2BED(norm_vcfs_input)
}



/*
// Process each VCF/CRAM pair separately
process subset_reads {
    publishDir "$params.outdir/mutect2_filtered_crams", mode: 'copy'
    input:
    path(vcf_file)
    path(cram_file)
    path(patient_id)
    path(sample_id)

    output:
    path("${patient_id}_${sample_id}.cram")

    script:
    """
    # Index the CRAM file
    samtools index ${cram_path}/${cram_file}

    # Extract the list of genomic positions from the VCF file
    bcftools query -f '%CHROM:%POS\n' ${vcf_path}/${vcf_file} > positions.txt

    # Subset the CRAM file to only include reads overlapping the positions in the VCF file
    samtools view -L positions.txt -b ${cram_path}/${cram_file} \
        | samtools sort -O CRAM -T tmp -o ${patient_id}_${sample_id}.cram \
        && samtools index ${patient_id}_${sample_id}.cram
    """
}
*/
