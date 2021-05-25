#!/bin/bash

conda create -n bwa -y bwa=0.7.17 samtools=1.12
conda activate bwa

mkdir -p data
cd data

# HG002
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/ILMN/downsampled/HG002_HiSeq30x_subsampled_R1.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/ILMN/downsampled/HG002_HiSeq30x_subsampled_R2.fastq.gz
# HG003
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/parents/ILMN/downsampled/HG003/HG003_HiSeq30x_subsampled_R1.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/parents/ILMN/downsampled/HG003/HG003_HiSeq30x_subsampled_R2.fastq.gz
# HG004
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/parents/ILMN/downsampled/HG004/HG004_HiSeq30x_subsampled_R1.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/parents/ILMN/downsampled/HG004/HG004_HiSeq30x_subsampled_R2.fastq.gz

mkdir -p ref
cd ref

# hg19 reference genome
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
gunzip human_g1k_v37.fasta.gz
bgzip human_g1k_v37.fasta
bwa index human_g1k_v37.fasta.gz

cd ..

bwa mem \
	-t 32 \
	-T 0 \
	-K 100000000 \
	-R '@RG\tID:HG002\tSM:HG002' \
	ref/human_g1k_v37.fasta.gz \
	HG002_HiSeq30x_subsampled_R1.fastq.gz \
	HG002_HiSeq30x_subsampled_R2.fastq.gz | \
	samtools sort --threads 32 --write-index -O BAM -o HG002.bam - 
samtools index HG002.bam
samtools view -b HG002.bam 22 X Y > HG002.22XY.bam

bwa mem \
	-t 32 \
	-T 0 \
	-K 100000000 \
	-R '@RG\tID:HG003\tSM:HG003' \
	ref/human_g1k_v37.fasta.gz \
	HG003_HiSeq30x_subsampled_R1.fastq.gz \
	HG003_HiSeq30x_subsampled_R2.fastq.gz | \
	samtools sort --threads 32 --write-index -O BAM -o HG003.bam - 
samtools index HG003.bam
samtools view -b HG003.bam 22 X Y > HG003.22XY.bam

bwa mem \
	-t 32 \
	-T 0 \
	-K 100000000 \
	-R '@RG\tID:HG004\tSM:HG004' \
	ref/human_g1k_v37.fasta.gz \
	HG004_HiSeq30x_subsampled_R1.fastq.gz \
	HG004_HiSeq30x_subsampled_R2.fastq.gz | \
	samtools sort --threads 32 --write-index -O BAM -o HG004.bam - 
samtools index HG004.bam
samtools view -b HG004.bam 22 X Y > HG004.22XY.bam

cd ..

conda deactivate
