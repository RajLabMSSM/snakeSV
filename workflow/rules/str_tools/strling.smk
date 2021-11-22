rule strling_index:
	input:
		REFERENCE_FASTA
	output:
		REFERENCE_FASTA + ".str"
	conda:
		SNAKEDIR + "envs/strling.yaml"
	log: 
		OUT_FOLDER + "/log/str_discovery/strling_index.log"	
	shell:
		"strling index {input} "

rule strling_extract:
	input:
		ref_fa = REFERENCE_FASTA,
		ref_index = REFERENCE_FASTA + ".str",
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		strling_bin = OUT_FOLDER + "/str_discovery/strling/{sample}/{sample}.bin"
	conda:
		SNAKEDIR + "envs/strling.yaml"
	shell:
		"strling extract -f {input.ref_fa} {input.bam} {output.strling_bin} "

# Get calls per sample
rule strling_call1:
	input:
		ref_fa = REFERENCE_FASTA,
		bam = OUT_FOLDER + "/input/{sample}.bam",
		strling_bin = OUT_FOLDER + "/str_discovery/strling/{sample}/{sample}.bin"
	output:
		strling_dir = directory(OUT_FOLDER + "/str_discovery/strling/{sample}/str-results"),
		strling_genotype = OUT_FOLDER + "/str_discovery/strling/{sample}/str-results/{sample}-genotype.txt"
	conda:
		SNAKEDIR + "envs/strling.yaml"
	shell:
		"strling call --output-prefix {output.strling_dir}/{wildcards.sample} -f {input.ref_fa} {input.bam} {input.strling_bin} "

# Multiple sample joint set
rule strling_merge:
	input:
		ref_fa = REFERENCE_FASTA,
		strling_bins = expand(OUT_FOLDER + "/str_discovery/strling/{sample}/{sample}.bin", sample=participant_id)
	output:
		strling_dir = directory(OUT_FOLDER + "/str_joint_calling/strling"),
		strling_bounds = OUT_FOLDER + "/str_joint_calling/strling/joint-bounds.txt"
	conda:
		SNAKEDIR + "envs/strling.yaml"
	shell:
		"strling merge --output-prefix {output.strling_dir}/joint -f {input.ref_fa} {input.strling_bins}"

# Calls from the joint set
rule strling_call2:
	input:
		ref_fa = REFERENCE_FASTA,
		bam = OUT_FOLDER + "/input/{sample}.bam",
		strling_bin = OUT_FOLDER + "/str_discovery/strling/{sample}/{sample}.bin",
		strling_bounds = OUT_FOLDER + "/str_joint_calling/strling/joint-bounds.txt"
	output:
		strling_dir = directory(OUT_FOLDER + "/str_joint_calling/{sample}/str-results"),
		strling_genotype = OUT_FOLDER + "/str_joint_calling/{sample}/str-results/{sample}-genotype.txt",
		strling_unplaced = OUT_FOLDER + "/str_joint_calling/{sample}/str-results/{sample}-unplaced.txt"
	conda:
		SNAKEDIR + "envs/strling.yaml"
	shell:
		"strling call --output-prefix {output.strling_dir}/{wildcards.sample} -b {input.strling_bounds} -f {input.ref_fa} {input.bam} {input.strling_bin} "

rule strling_outliers:
	input:
		strling_genotypes = expand(OUT_FOLDER + "/str_joint_calling/{sample}/str-results/{sample}-genotype.txt", sample=participant_id),
		strling_unplaced = expand(OUT_FOLDER + "/str_joint_calling/{sample}/str-results/{sample}-unplaced.txt", sample=participant_id),
	output:
		strling_genotype = OUT_FOLDER + "/str_joint_calling/strling_outliers/STRs.tsv"
	conda:
		SNAKEDIR + "envs/strling.yaml"
	shell:
		"strling-outliers.py --genotypes {input.strling_genotypes} --unplaced {input.strling_unplaced} "

