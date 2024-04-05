#!/bin/bash
set -e
# remove these files if found, if not found print warning message
if [ -f hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED ]; then
    rm hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED
else
    echo "Warning: hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED not found"
fi

if [ -f TIL10_signature.txt ]; then
    rm TIL10_signature.txt
else
    echo "Warning: TIL10_signature not found"
fi