#!/bin/bash

INPUT_DIR=/mnt/data/wes_data_copy/
for x in ${INPUT_DIR}*.fastq; do echo "$(basename -- $x)" && cat "$x" | paste - - - - | cut -f 2 | tr -d '\n' | wc -c; done > count_bases.txt
