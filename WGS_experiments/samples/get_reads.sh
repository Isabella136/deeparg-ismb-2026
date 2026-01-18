#!/bin/bash

while IFS= read -r biosample; do
    srr=`elink -db Biosample -id $biosample -target sra \
    | efetch -format docsum \
    | xtract -self -pattern DocumentSummary -if PAIRED -element Run@acc \
    | tr '\n' ' '`
    echo $srr;
    mkdir -p $biosample;
    cd $biosample;
    fasterq-dump $srr;
    cd ..;
done < ./../real_samples.txt