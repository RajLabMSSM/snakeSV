#!/bin/bash

curDir=${PWD}

cd ${curDir}/data || exit 1

mkdir -p ${curDir}/data/annotation
mkdir -p ${curDir}/data/annotation/custom
cd ${curDir}/data/annotation/custom

git clone https://github.com/nottalexi/brain-cell-type-peak-files.git

cd ${curDir}/data/annotation/custom/brain-cell-type-peak-files/H3K27ac

sed $'s/\r//' LHX2_optimal_peak.H3K27.bed | awk -v OFS='\t' '{gsub("chr","",$1); print $0 "\t" "Ast_H3K27ac"}' > ${curDir}/data/annotation/custom/astrocytes_H3K27ac.bed

sed $'s/\r//' NeuN_optimal_peak.H3K27.bed | awk -v OFS='\t' '{gsub("chr","",$1); print $0 "\t" "Neu_H3K27ac"}' > ${curDir}/data/annotation/custom/neurons_H3K27ac.bed

sed $'s/\r//' Olig2_optimal_peak.H3K27.bed | awk -v OFS='\t' '{gsub("chr","",$1); print $0 "\t" "Oli_H3K27ac"}' > ${curDir}/data/annotation/custom/oligodendrocytes_H3K27ac.bed

sed $'s/\r//' PU1_optimal_peak.H3K27.bed | awk -v OFS='\t' '{gsub("chr","",$1); print $0 "\t" "Mic_H3K27ac"}' > ${curDir}/data/annotation/custom/microglia_H3K27ac.bed

cd ${curDir}
