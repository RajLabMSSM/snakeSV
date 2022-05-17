#!/bin/bash

curDir=${PWD}

cd ${curDir}/data || exit 1

mkdir -p ${curDir}/data/annotation
mkdir -p ${curDir}/data/annotation/gtf
cd ${curDir}/data/annotation/gtf

wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh37_mapping/gencode.v38lift37.annotation.gtf.gz

conda activate snakesv_env
zcat gencode.v38lift37.annotation.gtf.gz | sed 's/^chr//g' | bgzip -c > gencode.v38lift37.annotation.nochr.gtf.gz

cd ${curDir}
