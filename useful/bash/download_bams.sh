#!/bin/bash
# array contains sample names
# replace lugh:// with lugh://path/to/alignment/folder
array=( 3532  4315 4340 4344 )
for i in "${array[@]}"
do
	start_time=$SECONDS
	echo "$i"
	scp lugh://"$i"/STAR/"$i".Aligned.sortedByCoord.out.bam .
	elapsed=$(( SECONDS - start_time ))
	eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"

done


