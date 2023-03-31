#!/bin/bash

INPUT_DIR=/mnt/data/wes_data_copy/
for x in ${INPUT_DIR}*.fastq; do echo "$(basename -- $x)" && cat "$x" | paste - - - - | cut -f 2 | tr -d '\n' | wc -c; done > count_bases.txt

#"$(basename -- $read1)"
#for x in ${INPUT_DIR}*.fastq; do cat "$x" | paste - - - - | cut -f 2 | tr -d '\n' | wc -c && echo $x; done
