#!/bin/bash
# need pyega3 and sendmail in path
# usage: bash pyega_download_all.sh -o path_to_outdir -c full_path_to_credentials_file.json

while getopts ":o::c:" flag;
do
    case $flag in
        o) outdir="$OPTARG";;
	c) credentials="$OPTARG";;
    esac
done
#use older version of pyega3
source activate pyega2


mkdir -p $outdir
#; cd $outdir

task(){
   mkdir -p ${2}/${1}; cd ${2}/${1}
   pwd
   echo "credentials: $3"
   pyega3 -cf $3 --connections 30 fetch $1
   echo "Subject: download script $1 finished" | sendmail k.ryan45@nuigalway.ie
}

# try not in parallel - in series - as the multiple connections might be messing up things
for dataset in EGAD00001006144 EGAD00001003808 EGAD00001004810 EGAD00001005744; do
    task "$dataset" "$outdir" "$credentials"
done
