#!/bin/bash
# need pyega3 and sendmail in path
# usage: bash pyega_download_all.sh -o path_to_outdir -c full_path_to_credentials_file.json
# TODO add your email address
# TODO create conda environment called pyega with pyega3 and sendmail installed
# TODO add your datasets
while getopts ":o::c:" flag;
do
    case $flag in
        o) outdir="$OPTARG";;
	c) credentials="$OPTARG";;
    esac
done


source activate pyega


mkdir -p $outdir

task(){
   mkdir -p ${2}/${1}; cd ${2}/${1}
   pwd
   echo "credentials: $3"
   pyega3 -cf $3 --connections 30 fetch $1
   echo "Subject: download script $1 finished" | sendmail youremail@domain.com
}

# add your datasets here
for dataset in EGAD1 EGAD2 EGAD3 EGAD4; do
    task "$dataset" "$outdir" "$credentials" &
done
