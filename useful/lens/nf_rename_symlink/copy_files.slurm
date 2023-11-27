#!/bin/bash
#SBATCH --time=16:00:00
#SBATCH --job-name="copy"
#SBATCH --output=copy_rename.out
#SBATCH --mail-user=k.ryan45@universityofgalway.ie
#SBATCH --mail-type=ALL
#SBATCH -p normal

start_time=$SECONDS

source activate nf-core

cd /data/kryan/misc/useful/lens/nf_rename

nextflow run main.nf -c nextflow.config

elapsed=$(( SECONDS - start_time ))
eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"

