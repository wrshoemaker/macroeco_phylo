#!/bin/bash



#while read sample; do
#  echo "$sample"
#  /Users/williamrshoemaker/sratoolkit.3.0.0-mac64/bin/fasterq-dump.3.0.0 -p --outdir /Users/williamrshoemaker/GitHub/strain_macroecology/data/poyet/poyet_16s/fastq/ ${sample}
#done </Users/williamrshoemaker/GitHub/macroeco_phylo/data/dalbello/SraAccList.txt



#/Users/williamrshoemaker/edirect  fasterq-dump

for filename in /Users/williamrshoemaker/SRA/sra/*.sra; do
    echo $filename
    /Users/williamrshoemaker/sratoolkit.3.0.0-mac64/bin/fasterq-dump.3.0.0 --split-files ${filename} --outdir /Users/williamrshoemaker/GitHub/macroeco_phylo/data/dalbello/fastq
done
