#!/bin/bash

curDir=${PWD}

cd ${curDir}/data || exit 1

mkdir -p ${curDir}/data/sv_panel
cd ${curDir}/data/sv_panel

# HG002
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v1.0/assemblies/20200628_HHU_assembly-results_CCS_v12/assemblies/phased/v12_NA24385_hpg_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v1.0/assemblies/20200628_HHU_assembly-results_CCS_v12/assemblies/phased/v12_NA24385_hpg_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta 

conda deactivate
conda create -n minimap2 -y minimap2=2.20 bcftools=1.12
conda activate minimap2

minimap2 --paf-no-hit -a -x asm5 --cs -r2k -t 32 ../ref/human_g1k_v37.fasta.gz v12_NA24385_hpg_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta > v12_NA24385_hpg_pbsq2-ccs_1000-pereg.h1-un.racon-p2.sam
minimap2 --paf-no-hit -a -x asm5 --cs -r2k -t 32 ../ref/human_g1k_v37.fasta.gz v12_NA24385_hpg_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta > v12_NA24385_hpg_pbsq2-ccs_1000-pereg.h2-un.racon-p2.sam

samtools sort -m4G -@4 -o v12_NA24385_hpg_pbsq2-ccs_1000-pereg.h1-un.racon-p2.sorted.bam v12_NA24385_hpg_pbsq2-ccs_1000-pereg.h1-un.racon-p2.sam
samtools sort -m4G -@4 -o v12_NA24385_hpg_pbsq2-ccs_1000-pereg.h2-un.racon-p2.sorted.bam v12_NA24385_hpg_pbsq2-ccs_1000-pereg.h2-un.racon-p2.sam

samtools index v12_NA24385_hpg_pbsq2-ccs_1000-pereg.h1-un.racon-p2.sorted.bam
samtools index v12_NA24385_hpg_pbsq2-ccs_1000-pereg.h2-un.racon-p2.sorted.bam

conda deactivate
conda create -n svim -y svim=1.4.2 bcftools=1.12
conda activate svim

svim-asm diploid hg002 v12_NA24385_hpg_pbsq2-ccs_1000-pereg.h1-un.racon-p2.sorted.bam v12_NA24385_hpg_pbsq2-ccs_1000-pereg.h2-un.racon-p2.sorted.bam human_g1k_v37.fasta.gz
bgzip -c hg002/variants.vcf > hg002/variants.vcf.gz
bcftools index hg002/variants.vcf.gz
bcftools view hg002/variants.vcf.gz --regions 22,X,Y | bgzip -c > hg002/variants.22XY.vcf.gz
bcftools index hg002/variants.22XY.vcf.gz
zcat hg002/variants.22XY.vcf.gz > hg002/variants.22XY.vcf

conda deactivate
cd ${curDir}
