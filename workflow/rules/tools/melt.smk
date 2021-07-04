#ln -f -s ${PWD}/resources/melt/me_refs/1KGP_Hg19 ${PWD}/resources/melt/me_refs/Hg37
#ln -f -s ${PWD}/resources/melt/add_bed_files/1KGP_Hg19 ${PWD}/resources/melt/add_bed_files/Hg37
#ln -f -s $PWD/resources/melt/add_bed_files/Hg37/hg19.genes.bed $PWD/resources/melt/add_bed_files/Hg37/hg37.genes.bed
#ln -f -s $PWD/resources/melt/add_bed_files/Hg37/hg19.FL.bed $PWD/resources/melt/add_bed_files/Hg37/hg37.FL.bed
#ln -f -s $PWD/resources/melt/add_bed_files/Hg37/hg19.FL.bed.source $PWD/resources/melt/add_bed_files/Hg37/hg37.FL.bed.source

rule melt_Preprocessing:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		disc = OUT_FOLDER + "/input/{sample}.bam.disc",
		disc_bai = OUT_FOLDER + "/input/{sample}.bam.disc.bai",
		fq = OUT_FOLDER + "/input/{sample}.bam.fq"
	params:
		melt_path = LIB_DIR + "melt"
	conda:
		SNAKEDIR + "envs/melt.yaml"
	shell:
		"java -Xmx2G -jar {params.melt_path}/MELT.jar Preprocess \
		-bamfile {input.bam} \
		-h {REFERENCE_FASTA}; "

rule melt_IndivAnalysis:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai",
		disc = OUT_FOLDER + "/input/{sample}.bam.disc",
		disc_bai = OUT_FOLDER + "/input/{sample}.bam.disc.bai",
		fq = OUT_FOLDER + "/input/{sample}.bam.fq"
	output:
		res_bam = OUT_FOLDER + "/sv_discovery/melt/{sample}/{MEI}/{sample}.{MEI}.aligned.final.sorted.bam",
		res_hum_breaks = OUT_FOLDER + "/sv_discovery/melt/{sample}/{MEI}/{sample}.{MEI}.hum_breaks.sorted.bam"
	params:
		folder = directory(OUT_FOLDER + "/sv_discovery/melt/{sample}/{MEI}"),
		melt_path = LIB_DIR + "melt",
		melt_ref = LIB_DIR + "melt/me_refs/Hg" + config["REF_BUILD"] + "/{MEI}_MELT.zip" 
	conda:
		SNAKEDIR + "envs/melt.yaml"
	shell:
		"java -Xmx6G -jar {params.melt_path}/MELT.jar IndivAnalysis \
  		-c 30 \
  		-h {REFERENCE_FASTA} \
  		-bamfile {input.bam} \
  		-t {params.melt_ref} \
  		-w {params.folder}; "
		"rm -rf {params.folder}/*tmp; "

rule melt_GroupAnalysis:
	input:
		res_bams = expand(OUT_FOLDER + "/sv_discovery/melt/{sample}/{MEI}/{sample}.{MEI}.aligned.final.sorted.bam", sample=participant_id, MEI=["ALU", "SVA", "LINE1"])
	output:
		pre_geno = OUT_FOLDER + "/sv_discovery/melt/Results/{MEI}/{MEI}.pre_geno.tsv"
	params:
		prevDir = directory(OUT_FOLDER + "/sv_discovery/melt/"),
		folder = directory(OUT_FOLDER + "/sv_discovery/melt/Results/{MEI}/"),
		melt_path = LIB_DIR + "melt",
		melt_bed = LIB_DIR + "melt/add_bed_files/Hg" + config["REF_BUILD"] + "/hg" + config["REF_BUILD"] + ".genes.bed",
		melt_ref = LIB_DIR + "melt/me_refs/Hg" + config["REF_BUILD"] + "/{MEI}_MELT.zip"
	conda:
		SNAKEDIR + "envs/melt.yaml"
	shell:
		"java -Xmx25G -jar {params.melt_path}/MELT.jar GroupAnalysis \
        -discoverydir {params.prevDir} \
        -w {params.folder} \
        -t {params.melt_ref} \
        -h {REFERENCE_FASTA} \
        -n {params.melt_bed}; "

rule melt_Genotype:
	input:
		pre_geno = expand(OUT_FOLDER + "/sv_discovery/melt/Results/{MEI}/{MEI}.pre_geno.tsv", MEI=["ALU", "SVA", "LINE1"])
	output:
		sample_geno = OUT_FOLDER + "/sv_discovery/melt/{sample}/{MEI}.gt/{sample}.{MEI}.tsv"
	params:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		step2_folder = directory(OUT_FOLDER + "/sv_discovery/melt/Results/{MEI}/"),
		folder = directory(OUT_FOLDER + "/sv_discovery/melt/{sample}/{MEI}.gt/"),
		melt_path = LIB_DIR + "melt",
		melt_bed = LIB_DIR + "melt/add_bed_files/Hg" + config["REF_BUILD"] + "/hg" + config["REF_BUILD"] + ".genes.bed",
		melt_ref = LIB_DIR + "melt/me_refs/Hg" + config["REF_BUILD"] + "/{MEI}_MELT.zip"
	conda:
		SNAKEDIR + "envs/melt.yaml"
	shell:
		"java -Xmx2G -jar {params.melt_path}/MELT.jar Genotype \
  		-bamfile {params.bam} \
  		-t {params.melt_ref} \
  		-h {REFERENCE_FASTA} \
  		-w {params.folder} \
  		-p {params.step2_folder}; "
  		"rm -rf {params.folder}/*tmp; "

rule melt_MakeVCF:
	input:
		sample_geno = OUT_FOLDER + "/sv_discovery/melt/{sample}/{MEI}.gt/{sample}.{MEI}.tsv"
	output:
		vcf = OUT_FOLDER + "/sv_discovery/melt/{sample}/{MEI}.final_comp.vcf",
		vcfgz = OUT_FOLDER + "/sv_discovery/melt/{sample}/{MEI}.final_comp.vcf.gz",
		tbi = OUT_FOLDER + "/sv_discovery/melt/{sample}/{MEI}.final_comp.vcf.gz.tbi"
	params:
		root_folder = directory(OUT_FOLDER + "/sv_discovery/melt/{sample}/"),
		step2_folder = directory(OUT_FOLDER + "/sv_discovery/melt/Results/{MEI}/"),
		step3_folder = directory(OUT_FOLDER + "/sv_discovery/melt/{sample}/{MEI}.gt/"),
		melt_path = LIB_DIR + "melt",
		melt_ref = LIB_DIR + "melt/me_refs/Hg" + config["REF_BUILD"] + "/{MEI}_MELT.zip"
	conda:
		SNAKEDIR + "envs/melt.yaml"
	shell:
		"java -Xmx2G -jar {params.melt_path}/MELT.jar MakeVCF \
    	-genotypingdir {params.step3_folder} \
    	-h {REFERENCE_FASTA} \
    	-t {params.melt_ref} \
    	-w {params.root_folder} \
    	-p {params.step2_folder} \
    	-o {params.root_folder}; "
    	"bgzip -c {output.vcf} > {output.vcfgz}; "
    	"tabix -p vcf {output.vcfgz}; "

rule melt_MergeVCF:
	input:
		vcfs = lambda wildcards: expand(OUT_FOLDER + "/sv_discovery/melt/" + wildcards.sample + "/{MEI}.final_comp.vcf", MEI=["ALU", "SVA", "LINE1"])
	output:
		vcf = OUT_FOLDER + "/sv_discovery/melt/{sample}/{sample}.melt.vcf.gz",
		tbi = OUT_FOLDER + "/sv_discovery/melt/{sample}/{sample}.melt.vcf.gz.tbi"
	conda:
		SNAKEDIR + "envs/melt.yaml"
	shell:
		"bcftools concat -a $(find {input.vcfs} -maxdepth 1 -size +0 -print | awk '{{print $1\".gz\"}}') | bcftools sort -Oz -o {output.vcf}; "
		"tabix -p vcf {output.vcf}; "