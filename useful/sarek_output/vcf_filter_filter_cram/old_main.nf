#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
// Read in the list of VCF files
vcf_files = Channel
    .fromPath(params.vcf_path)
    .map { file -> file.getName() }
    //.view()
// Read in the list of CRAM files
cram_files = Channel
    .fromPath(params.cram_path)
    .map { file -> file.getName() }
    //.view()
// Combine the VCF and CRAM files based on the sample IDs
vcf_cram_pairs = vcf_files
    .combine(cram_files)
    .view()
    .map { pair ->
        def sample_ids = pair[0].tokenize('_vs_')
        def cram_file = cram_files.find { file -> file.startsWith(sample_ids[0]) }
        [pair[0], cram_file, sample_ids[0], sample_ids[1]]
    }
    .view()
*/

/*
vcf_files = Channel.fromPath(params.vcf_path)
                   //.view()

// Read in the list of CRAM files
cram_files = Channel.fromPath(params.cram_path)
                  //  .view()
// Combine the VCF and CRAM files based on the sample IDs
/*
vcf_cram_pairs = vcf_files
    .mix(cram_files)
    .map { pair ->
        def sample_ids = pair[0].tokenize('_vs_')
        def cram_file = cram_files.find { file -> file.getName().startsWith(sample_ids[0]) }
        [pair[0], cram_file.getName(), sample_ids[0], sample_ids[1]]
    }
*/
/*
vcf_cram_pairs = vcf_files
    .mix(cram_files)
    .map { pair ->
        def sample_ids = pair[0].toString().tokenize('_vs_')
        def cram_file = cram_files.find { file -> file.getName().startsWith(sample_ids[0]) }
        [pair[0], cram_file.getName(), sample_ids[0], sample_ids[1]]
    }
    .view()
*/

// Read in the list of VCF files
vcf_files = Channel.fromPath(params.vcf_path)
    .map { file(it).getName() } // Convert to file objects and get the name
    //.view()
// Read in the list of CRAM files
cram_files = Channel.fromPath(params.cram_path)
    .map { file(it).getName() } // Convert to file objects and get the name
    //.view()
// Combine the VCF and CRAM files based on the sample IDs
vcf_cram_pairs = vcf_files
    .mix(cram_files)
    .map { pair ->
        def sample_ids = pair[0].getName().tokenize('_vs_')
        def cram_file = cram_files.find { file(it).getName().startsWith(sample_ids[0]) }
        [pair[0], cram_file.getName(), sample_ids[0], sample_ids[1]]
    }

vcf_cram_pairs.view()
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
