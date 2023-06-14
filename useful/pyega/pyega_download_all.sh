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

source activate pyega


mkdir -p $outdir; cd $outdir
mkdir -p EGAD00001006144; cd EGAD00001006144
pyega3 -cf $credentials fetch EGAD00001006144

echo "Subject: download script EGAD00001006144 finished" | sendmail k.ryan45@nuigalway.ie

cd ../
mkdir -p EGAD00001003808; cd EGAD00001003808
pyega3 -cf $credentials fetch EGAD00001003808

echo "Subject: download script EGAD00001003808 finished" | sendmail k.ryan45@nuigalway.ie

# size: 143 Gb
cd ../
mkdir -p EGAD00001004810; cd EGAD00001004810
pyega3 -cf $credentials fetch EGAD00001004810
echo "Subject: download script EGAD00001004810 finished" | sendmail k.ryan45@nuigalway.ie

cd ../
mkdir -p EGAD00001005744; cd EGAD00001005744
# size: 94Gb
pyega3 -cf $credentials fetch EGAD00001005744
echo "Subject: download script EGAD00001005744 finished" | sendmail k.ryan45@nuigalway.ie

