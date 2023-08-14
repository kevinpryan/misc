#!/bin/bash
#SBATCH --time=16:00:00
#SBATCH --job-name="copy"
#SBATCH --output=copy_rename.out
#SBATCH --mail-user=<email>
#SBATCH --mail-type=ALL
#SBATCH -p normal

start_time=$SECONDS

source activate nf-core

cd /data2/kryan/lens/

nextflow run main.nf -c nextflow.config

elapsed=$(( SECONDS - start_time ))
eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"

