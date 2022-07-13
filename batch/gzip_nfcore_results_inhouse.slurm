#!/bin/bash
#SBATCH --time=16:00:00
#SBATCH --job-name="gzip"
#SBATCH --output=gzip.out
#SBATCH --mail-user=your_email_address
#SBATCH --mail-type=ALL
#SBATCH -N 1 # number of nodes

start_time=$SECONDS

cd /data2/kryan/inhouse_rnaseq/
# usage: tar -zcvf name_file.tar.gz output_file_name/
tar -zcvf inhouse_data_nfcore_results.tar.gz inhouse_data_nfcore_results/

elapsed=$(( SECONDS - start_time ))
eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"
