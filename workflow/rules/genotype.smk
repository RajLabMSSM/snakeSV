###################################################################################################
## Merge tools
###################################################################################################
rule decompress_for_merging:
	input:
		vcf = OUT_FOLDER + "/sv_discovery/{tool}/{sample}/{sample}.{tool}.vcf.gz"
	output:
		vcf = OUT_FOLDER + "/sv_discovery/merged/{sample}/{sample}.{tool}.vcf"
	shell:	
		"zcat {input.vcf} > {output.vcf}"

rule merge_tools:
	input:
		vcf = expand(OUT_FOLDER + "/sv_discovery/merged/{{sample}}/{{sample}}.{tool}.vcf", tool=TOOLS)
	output:
		vcf = OUT_FOLDER + "/sv_discovery/merged/{sample}/{sample}.vcf",
		vcfgz = OUT_FOLDER + "/sv_discovery/merged/{sample}/{sample}.vcf.gz",
		tbi = OUT_FOLDER + "/sv_discovery/merged/{sample}/{sample}.vcf.gz.tbi"
	conda:
		"../envs/survivor.yaml"
	params:
		vcf_list = OUT_FOLDER + "/sv_discovery/merged/{sample}/vcf.list",
		breakpoint_dist = "500",
		min_num_calls = "1",
		use_type = "1",
		use_strand = "1",
		dist_based = "0",
		min_sv_size = "50"
	shell:
		"ls {input.vcf} > {params.vcf_list}; "
		#"svimmer {params.vcf_list} \{1..22\} > {output.vcf}; "
		"SURVIVOR merge {params.vcf_list} \
			{params.breakpoint_dist} \
			{params.min_num_calls} \
			{params.use_type} \
			{params.use_strand} \
			{params.dist_based} \
			{params.min_sv_size} \
			{output.vcf}; "
		"bcftools sort -Oz -o {output.vcfgz} {output.vcf}; "
		"tabix -p vcf {output.vcfgz}; "

###################################################################################################
## Joint-genotyping
###################################################################################################
rule merge_samples_01:
	input:
		vcf = expand(OUT_FOLDER + "/sv_discovery/merged/{sample}/{sample}.vcf", sample=participant_id)
	output:
		vcf_list = OUT_FOLDER + "/merged_cohort/vcf.list"
	run:
		shell("ls {input.vcf} > {output.vcf_list}; ")
		if( "SV_PANEL" in config ):
				for sv_ref in config["SV_PANEL"]:
					shell("ls {sv_ref} >> {output.vcf_list}; ")

rule merge_samples_02:
	input:
		vcf_list = OUT_FOLDER + "/merged_cohort/vcf.list"
	output:
		vcf = OUT_FOLDER + "/merged_cohort/raw_merged.vcf",
		vcfgz = OUT_FOLDER + "/merged_cohort/raw_merged.vcf.gz",
		tbi = OUT_FOLDER + "/merged_cohort/raw_merged.vcf.gz.tbi"
	conda:
		"../envs/survivor.yaml"
	params:
		vcf_list = OUT_FOLDER + "/merged_cohort/vcf.list",
		breakpoint_dist = "500",
		min_num_calls = "1",
		use_type = "1",
		use_strand = "1",
		dist_based = "0",
		min_sv_size = "50"
	shell:
		"SURVIVOR merge {params.vcf_list} \
			{params.breakpoint_dist} \
			{params.min_num_calls} \
			{params.use_type} \
			{params.use_strand} \
			{params.dist_based} \
			{params.min_sv_size} \
			{output.vcf}; "
		"bcftools sort -Oz -o {output.vcfgz} {output.vcf}; "
		"tabix -p vcf {output.vcfgz}; "

rule genotype:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		gvcf = OUT_FOLDER + "/merged_cohort/raw_merged.vcf.gz"
	output:
		vcf = OUT_FOLDER + "/sv_genotyping/{sample}/{sample}.vcf.gz",
		tbi = OUT_FOLDER + "/sv_genotyping/{sample}/{sample}.vcf.gz.tbi"
	conda:
		"../envs/graphtyper.yaml"
	params:
		outdir = directory(OUT_FOLDER + "/sv_genotyping/{sample}")
	shell:
		"for chr in $(tabix -l {input.gvcf}); do \
			graphtyper genotype_sv {REFERENCE_FASTA} {input.gvcf} \
			--log={params.outdir}/$chr/{wildcards.sample}.log \
			--sam={input.bam} \
			--region=$chr \
			--output={params.outdir}; \
		done; "
		"graphtyper vcf_concatenate --no_sort {params.outdir}/*/*.vcf.gz | bgzip -c > {output.vcf};"
		"tabix -p vcf {output.vcf}"

rule merge_genotyped_samples:
	input:
		vcf = expand(OUT_FOLDER + "/sv_genotyping/{sample}/{sample}.vcf.gz", sample=participant_id)
	output:
		vcf = OUT_FOLDER + "/merged_cohort/gt_merged.vcf",
		vcfgz = OUT_FOLDER + "/merged_cohort/gt_merged.vcf.gz",
		tbi = OUT_FOLDER + "/merged_cohort/gt_merged.vcf.gz.tbi"
	conda:
		"../envs/graphtyper.yaml"		
	shell:
		"graphtyper vcf_merge {input.vcf} --sv | grep -E -v \"BREAKPOINT|COVERAGE\" > {output.vcf}; "
		"bcftools sort -Oz -o {output.vcfgz} {output.vcf}; "
		"tabix -p vcf {output.vcfgz}; "
