#!/bin/bash

echo "trying after adding the iedb_install_ok file, running on full files"

NXF_VER=20.10.0 nextflow -bg run nextNEOpi.nf \
    --batchFile batchFile.csv \
    -config conf/params.config \
    --outputDir nextNeoPi_results \
    --trim_adapters false \
    --trim_adapters_RNAseq false \
    --use_NetChop false \
    --accept_license \
    -N k.ryan45@nuigalway.ie \
    -profile singularity

