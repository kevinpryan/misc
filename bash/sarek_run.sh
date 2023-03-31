nextflow run nf-core-sarek-3.1.2/workflow/main.nf -profile singularity --input /home/administrator/Documents/wes_samplesheet.csv --outdir /mnt/data/sarek_out/ --wes --intervals /mnt/data/references/Twist_Comprehensive_Exome_Covered_Targets_hg38_sorted.bed --only_paired_variant_calling --genome GATK.GRCh38 --email k.ryan45@nuigalway.ie --save_reference --max_cpus 20 --max_memory 125.5GB --tools strelka,mutect2,vep,snpeff -resume

nextflow run nf-core-sarek-3.1.2/workflow/main.nf -profile singularity --step variant_calling --input /mnt/data/sarek_out/csv/fixed_recalibrated.csv --outdir /mnt/data/variant_calling_sarek_out/ --wes --intervals /mnt/data/references/Twist_Comprehensive_Exome_Covered_Targets_hg38_sorted.bed --only_paired_variant_calling --genome GATK.GRCh38 --email k.ryan45@nuigalway.ie --save_reference --max_cpus 20 --max_memory 125.5GB --tools strelka,mutect2,freebayes,manta,tiddit,msisensorpro,cnvkit,controlfreec,merge -resume

# germline variant calling
nextflow run nf-core-sarek-3.1.2/workflow/main.nf -profile singularity --step variant_calling --input /mnt/data/sarek_out/csv/fixed_recalibrated.csv --outdir /mnt/data/germline_variant_calling_sarek_out/ --wes --intervals /mnt/data/references/Twist_Comprehensive_Exome_Covered_Targets_hg38_sorted.bed --genome GATK.GRCh38 --email k.ryan45@nuigalway.ie --save_reference --max_cpus 20 --max_memory 125.5GB --tools strelka,freebayes,deepvariant,haplotypecaller,mpileup,merge --joint_germline -resume

