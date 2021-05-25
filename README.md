# *snakeSV*: Flexible framework for large-scale SV discovery

Ricardo A. Vialle (ricardovialle@gmail.com)

---

*snakeSV* is an integrated pipeline in Snakemake for complete SV analysis. The pipeline includes pre- and post-processing steps to deal with large scale studies. The input data of the pipeline consists of BAM files for each sample, a reference genome file (.FASTA) and a configuration file in yaml format. Additionally, users can also input custom annotation files in BED format for SV interpretation and VCF files with structural variants to be genotyped in addition to the discovery set.

![Pipeline Schematic](docs/Pipeline_Schema.png "Pipeline Schematic")

---

### Requirements:

The mandatory requirements are a Linux environment with Python and Git. The pipeline focuses on using Conda environments to centralize all the tool management, but it can be easily customizable to include different tools and methods, not necessarily distributed by Anaconda (and derivatives).

#### 1-Miniconda:

This step can be ignored if Anaconda is already installed in the system.
```
# Download conda installer
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Set permissions to execute
chmod +x Miniconda3-latest-Linux-x86_64.sh 	

# Execute. Make sure to "yes" to add the conda to your PATH
./Miniconda3-latest-Linux-x86_64.sh 		

# Add channels
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```

#### 2-Snakemake:

Install Snakemake:
```
conda install snakemake
```

Alternatively, install Snakemake in a separated environment (named "snakemake_env") with the command:
```
conda create -n snakemake_env snakemake
conda activate snakemake_env # Command to activate the environment. To deactivate use "conda deactivate"
```

#### 3-snakeSV:

To clone the repository and install the required tools and download required supporting data use the following commands:
```
git clone https://github.com/RajLabMSSM/snakeSV.git
cd snakeSV
```

#### 4-Resources:

This repository comes with files required for running Delly, Smoove and somalier. 

* Delly blacklist regions: [downloaded from here](https://gear.embl.de/data/delly/)
* Smoove regions to exclude (from SpeedSeq): 
[GRCh37](https://github.com/hall-lab/speedseq/blob/master/annotations/ceph18.b37.lumpy.exclude.2014-01-15.bed) |
[GRCh38]( https://github.com/hall-lab/speedseq/blob/master/annotations/exclude.cnvnator_100bp.GRCh38.20170403.bed)
* Somalier sites: 
[sites.hg19.vcf.gz](https://github.com/brentp/somalier/files/3412453/sites.hg19.vcf.gz) |
[sites.hg38.nochr.vcf.gz](https://github.com/brentp/somalier/files/3412454/sites.hg38.nochr.vcf.gz) |
[sites.GRCh37.vcf.gz](https://github.com/brentp/somalier/files/3412455/sites.GRCh37.vcf.gz) |
[sites.hg38.vcf.gz](https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz)

Decompress it in snakeSV folder (a folder named resources should be created)

```
tar -zxvf resources.tar.gz
```

#### 5-Tools:

All other tools and dependencies are installed in their own environment automatically on the first run (with `--use-conda` parameter active). 

```
snakemake --configfile config/config.yaml --use-conda --conda-create-envs-only --cores 2
```

---

### Run:

#### Example data:

An example test run can be performed using the whole-genome data from the Ashkenazi Jewish trio (HG002-son, HG003-father, and HG004-mother) from the The Personal Genome Project. Sequencing for this data was performed using the Illumina platform with PCR-free and 150 bp pair-ended and data was downsampled to ~30x to match the expected production from most available studies (e.g. 1000 genome samples). Links for downloading raw files for each sample are available via the Human Pangenome Reference Consortium (https://github.com/human-pangenomics/HG002_Data_Freeze_v1.0) in the following links: 

* [HG002](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/ILMN/downsampled/)
* [HG003 and HG004](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/parents/ILMN/downsampled/)

Download and alignment of these samples can be performed using the following script (files will saved at `data` folder):
```
example/01_prepare_short_read.sh
```

SVs discovered using long reads, can also be included for genotyping using short-reads. Here we will use diploid assemblies for the HG002 sample, generated from the Human Genome Structural Variation Consortium (HGSVC) (Ebert et al. 2021). Links to each haplotype assemblies [H1](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v1.0/assemblies/20200628_HHU_assembly-results_CCS_v12/assemblies/phased/v12_NA24385_hpg_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta) adnd [H2](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v1.0/assemblies/20200628_HHU_assembly-results_CCS_v12/assemblies/phased/v12_NA24385_hpg_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta)

Similarly as before, assembly sequences, alignment and SV discovery can be performed using the following script:
```
example/02_prepare_long_read.sh
```

Cell-specific brain epigenetic annotations generated from ChIP-seq and ATAC-seq experiments (Nott et al. 2019). Data can be downloaded from the following GitHub repository:
https://github.com/nottalexi/brain-cell-type-peak-files

Script for downloading these annotations, and adapting files for snakeSV run:
```
example/03_prepare_annotations.sh
```

#### ConfigFile:

To perform an analysis, the user needs to specify a field SAMPLE_KEY in the yaml configuration pointing to a file with a unique identification for each sample and the paths to the respective BAM files. The SAMPLE_KEY file must be tabulated containing two columns named “participant_id” and “bam”.  sampleKey file for the example is found at `example/sampleKey.txt`

```
participant_id	bam
HG002	data/HG002.22XY.bam
HG003	data/HG003.22XY.bam
HG004	data/HG004.22XY.bam
```

Additionally, fields pointing to the reference genome fasta file (REFERENCE_FASTA), the reference build used (37 or 38, REF_BUILD), the output folder to store the results (OUT_FOLDER), and the list of SV discovery tools to be used (TOOLS) must be specified. Other supplementary files are specified by default with data distributed in the repository, but can also be customized. Fields ANNOTATION_BED and SV_PANEL are optional, and will be explained better in the use case section. Note that all these parameters can also be specified through the command-line interface, using the snakemake -C argument (-C [KEY=VALUE [KEY=VALUE ...]]). An example of a configuration file is shown below, or [here](config/config.yaml).

```
OUT_FOLDER: "results1"
SAMPLE_KEY: "example/sampleKey.txt"

TOOLS: 
  - "manta"
  - "smoove" 
  - "delly"

ANNOTATION_BED: 
  - "data/brain-cell-type-peak-files/astrocytes_H3K27ac.bed"
  - "data/brain-cell-type-peak-files/microglia_H3K27ac.bed"
  - "data/brain-cell-type-peak-files/neurons_H3K27ac.bed"
  - "data/brain-cell-type-peak-files/oligodendrocytes_H3K27ac.bed"

SV_PANEL: 
  - "data/sv_panel/hg002/variants.22.vcf"

LIB_DIR: "resource/"

REFERENCE_FASTA: "data/ref/human_g1k_v37.fasta"
REF_BUILD: "37"
GENCODE_GTF: "data/annotation/gencode.v38lift37.annotation.nochr.gtf.gz"
DICT: "data/ref/human_g1k_v37.dict"
NMASK: "data/ref/human_g1k_v37.mask.36.fasta.bed"

TMP_DIR: "~/ad-omics/ricardo/tmp/"
```

After having all the required data done, perform the analysis as follows:
```
snakemake --configfile config/config.yaml -pr --cores 2 --use-conda
```

---

### HPC run:

Make a copy of cluster configuration file:
```
cp config/cluster_lsf.yaml cluster.yaml
```

Edit the file with your cluster specifications (threads, partitions, cpu/memory, etc) for each rule.

Run snakeSV (LSF example):
```
./snakejob -u cluster.yaml -c config/config.yaml
```
