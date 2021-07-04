###################################################################################################
## Benchmarking against GiaB
###################################################################################################
# Deactivated in the pipeline.
rule benchmark_raw:
	input:
		vcf = OUT_FOLDER + "/sv_discovery/{tool}/{sample}/{sample}.{tool}.vcf.gz"
	output:
		OUT_FOLDER + "/bench/{sample}/sv_discovery/{tool}/giab_report.txt"
	conda:
		SNAKEDIR + "envs/truvari.yaml"		
	params:
		base = "data/bench/HG002_SVs_Tier1_v0.6.vcf.gz",
		bed = "data/bench/HG002_SVs_Tier1_v0.6.22XY.bed",
		ref = REFERENCE_FASTA,
		outdir = OUT_FOLDER + "/bench/{sample}/sv_discovery/{tool}/"
	shell:
		"rm -rf {params.outdir}; "
		"truvari bench \
 			-b {params.base} \
 			-c {input.vcf} \
 			-o {params.outdir} \
 			-f {params.ref} \
 			--includebed {params.bed} \
 			--passonly \
 			-p 0 \
 			-r 2000 \
 			--giabreport; "

rule benchmark_raw_gt:
	input:
		vcf = OUT_FOLDER + "/sv_discovery/{tool}/{sample}/{sample}.{tool}.gt.vcf.gz"
	output:
		OUT_FOLDER + "/bench/{sample}/sv_discovery/{tool}.gt/giab_report.txt"
	params:
		base = "data/bench/HG002_SVs_Tier1_v0.6.vcf.gz",
		bed = "data/bench/HG002_SVs_Tier1_v0.6.22XY.bed",
		ref = REFERENCE_FASTA,
		outdir = OUT_FOLDER + "/bench/{sample}/sv_discovery/{tool}.gt/"
	conda:
		SNAKEDIR + "envs/truvari.yaml"
	shell:
		"rm -rf {params.outdir}; "
		"zcat {input.vcf} | grep -E \"^#|AGGREGATED\" | bcftools sort -Oz -o {input.vcf}.filt.vcf.gz; "
		"tabix -p vcf {input.vcf}.filt.vcf.gz; "
		"truvari bench \
 			-b {params.base} \
 			-c {input.vcf}.filt.vcf.gz \
 			-o {params.outdir} \
 			-f {params.ref} \
 			--includebed {params.bed} \
 			--passonly \
 			-p 0 \
 			-r 2000 \
 			--giabreport; "

rule benchmark_gt:
	input:
		vcf = OUT_FOLDER + "/sv_genotyping/{sample}/{sample}.vcf.gz"
	output:
		OUT_FOLDER + "/bench/{sample}/sv_genotyping/{tool}/giab_report.txt"
	params:
		base = "data/bench/HG002_SVs_Tier1_v0.6.vcf.gz",
		bed = "data/bench/HG002_SVs_Tier1_v0.6.22XY.bed",
		ref = REFERENCE_FASTA,
		outdir = OUT_FOLDER + "/bench/{sample}/sv_genotyping/{tool}/"
	conda:
		SNAKEDIR + "envs/truvari.yaml"
	shell:
		"rm -rf {params.outdir}; "
		"zcat {input.vcf} | grep -E \"^#|AGGREGATED\" | bcftools sort -Oz -o {input.vcf}.filt.vcf.gz; "
		"tabix -p vcf {input.vcf}.filt.vcf.gz; "
		"truvari bench \
 			-b {params.base} \
 			-c {input.vcf}.filt.vcf.gz \
 			-o {params.outdir} \
 			-f {params.ref} \
 			--includebed {params.bed} \
 			--passonly \
 			-p 0 \
 			-r 2000 \
 			--giabreport; "

rule benchmark_merged_gt:
	input:
		vcf = OUT_FOLDER + "/merged_cohort/gt_merged.vcf.gz"
	output:
		OUT_FOLDER + "/bench/{sample}/merged/giab_report.txt"
	params:
		base = "data/bench/HG002_SVs_Tier1_v0.6.vcf.gz",
		bed = "data/bench/HG002_SVs_Tier1_v0.6.22XY.bed",
		ref = REFERENCE_FASTA,
		outdir = OUT_FOLDER + "/bench/{sample}/merged/"
	conda:
		SNAKEDIR + "envs/truvari.yaml"
	shell:
		"rm -rf {params.outdir}; "
		"truvari bench \
 			-b {params.base} \
 			-c {input.vcf} \
 			-o {params.outdir} \
 			-f {params.ref} \
 			--includebed {params.bed} \
 			-p 0 \
 			-r 2000 \
 			--giabreport; "
