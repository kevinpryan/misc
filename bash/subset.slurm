#!/bin/bash
#SBATCH --time=16:00:00
#SBATCH --job-name="subset"
#SBATCH --output=subset_rna.out
#SBATCH --mail-user=<email>
#SBATCH --mail-type=ALL
#SBATCH -p normal

start_time=$SECONDS

# subset reads giving a csv as input, cols old, new representing filenames. seqtk can subsample gzipped files, but output is not gzipped - have to gzip them yourself
# usage: subset.sh -f sub.txt -d $input_dir -o $output_dir

while getopts f:d:o: flag
do
    case "${flag}" in
        f) new_filenames=${OPTARG};;
        d) input_dir=${OPTARG};;
        o) output_dir=${OPTARG};;
    esac
done

cd $input_dir

    while IFS="," read -r old new; do
	echo $old
        echo $new
        OUT=${output_dir}/${new}
        singularity run -B <dir-to-bind> <path-to-seqtk-image> seqtk sample -s100 $old 10000 > $OUT
        gzip -f $OUT
    done < <(tail -n +2 $new_filenames)

elapsed=$(( SECONDS - start_time ))
eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"

