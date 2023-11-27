#!/bin/bash

# use samplesheet output of rename_files_general.R
# usage: rename.sh -f samplesheet.csv -d $DNA_DIR -r $RNA_DIR

while getopts f:d:r: flag
do
    case "${flag}" in
        f) new_filenames=${OPTARG};;
        d) input_dir_dna=${OPTARG};;
        r) input_dir_rna=${OPTARG};;
    esac
done

    while IFS="," read -r file_original patient abnormal_normal file_lens Sequencing_Method; do
        # dna and rna samples are in different directories. 
        if [ $Sequencing_Method == 'WES' ] || [ $Sequencing_Method == 'WGS' ] || [ $Sequencing_Method == 'WXS' ]
        then
            ORIGINAL_PATH="${input_dir_dna}/${file_original}"
        elif [ $Sequencing_Method == "RNA" ]
        then
            ORIGINAL_PATH="${input_dir_rna}/${file_original}"
        fi
        if [ -f $ORIGINAL_PATH ] ; then
                echo "original path: $ORIGINAL_PATH"
                echo "new name: $file_lens"
                #copying and renaming at the same time
                cp $ORIGINAL_PATH $file_lens
                #mv $NR $NAME
        fi
    done < <(tail -n +2 $new_filenames)
    
