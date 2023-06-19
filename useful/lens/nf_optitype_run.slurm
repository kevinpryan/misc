#!/bin/bash
#SBATCH --time=10-00:00
#SBATCH --job-name="optitype"
#SBATCH --output=optitype.out
#SBATCH --mail-user=k.ryan45@nuigalway.ie
#SBATCH --mail-type=ALL
#SBATCH -p normal

start_time=$SECONDS

source activate nf-core
OPTITPE_PATH="/data4/kryan/misc/useful/lens/nf_optitype"
cd $OPTITYPE_PATH
nextflow run main.nf -c nextflow.config                                                                                                 

elapsed=$(( SECONDS - start_time ))
eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"

