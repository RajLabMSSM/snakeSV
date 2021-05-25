#!/bin/bash

cd data

mkdir -p annotation
cd annotation

wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/GRCh37_mapping/gencode.v38lift37.annotation.gtf.gz
zcat gencode.v38lift37.basic.annotation.gtf.gz | sed 's/^chr//g' | bgzip -c > gencode.v38lift37.annotation.nochr.gtf.gz

cd ..

git clone https://github.com/nottalexi/brain-cell-type-peak-files.git

cd brain-cell-type-peak-files/H3K27ac

sed $'s/\r//' LHX2_optimal_peak.H3K27.bed | awk -v OFS='\t' '{gsub("chr","",$1); print $0 "\t" "Ast_H3K27ac"}' > ../astrocytes_H3K27ac.bed

sed $'s/\r//' NeuN_optimal_peak.H3K27.bed | awk -v OFS='\t' '{gsub("chr","",$1); print $0 "\t" "Neu_H3K27ac"}' > ../neurons_H3K27ac.bed

sed $'s/\r//' Olig2_optimal_peak.H3K27.bed | awk -v OFS='\t' '{gsub("chr","",$1); print $0 "\t" "Oli_H3K27ac"}' > ../oligodendrocytes_H3K27ac.bed

sed $'s/\r//' PU1_optimal_peak.H3K27.bed | awk -v OFS='\t' '{gsub("chr","",$1); print $0 "\t" "Mic_H3K27ac"}' > ../microglia_H3K27ac.bed

cd ../../..

