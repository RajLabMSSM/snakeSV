###################################################################################################
## Benchmarking
###################################################################################################
rule benchmark_raw:
	input:
		vcf = OUT_FOLDER + "/sv_discovery/{tool}/{sample}/{sample}.{tool}.vcf.gz"
	output:
		OUT_FOLDER + "/sv_discovery/{tool}/{sample}/bench/giab_report.txt"
	conda:
		"../envs/truvari.yaml"		
	params:
		base = "~/ad-omics/ricardo/Data/GIAB/HG002_SV_benchmark/HG002_SVs_Tier1_v0.6.vcf.gz",
		bed = "~/ad-omics/ricardo/Data/GIAB/HG002_SV_benchmark/HG002_SVs_Tier1_v0.6.bed",
		ref = REFERENCE_FASTA,
		outdir = OUT_FOLDER + "/sv_discovery/{tool}/{sample}/bench/"
	conda:
		"../envs/truvari.yaml"
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
		OUT_FOLDER + "/sv_discovery/{tool}/{sample}/bench_gt/giab_report.txt"
	params:
		base = "~/ad-omics/ricardo/Data/GIAB/HG002_SV_benchmark/HG002_SVs_Tier1_v0.6.vcf.gz",
		bed = "~/ad-omics/ricardo/Data/GIAB/HG002_SV_benchmark/HG002_SVs_Tier1_v0.6.bed",
		ref = REFERENCE_FASTA,
		outdir = OUT_FOLDER + "/sv_discovery/{tool}/{sample}/bench_gt/"
	conda:
		"../envs/truvari.yaml"
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
 			-p 0 \
 			-r 2000 \
 			--giabreport; "


rule benchmark_gt:
	input:
		vcf = OUT_FOLDER + "/sv_genotyping/{sample}/{sample}.vcf.gz"
	output:
		OUT_FOLDER + "/sv_genotyping/{sample}/bench/giab_report.txt"
	params:
		base = "~/ad-omics/ricardo/Data/GIAB/HG002_SV_benchmark/HG002_SVs_Tier1_v0.6.vcf.gz",
		bed = "~/ad-omics/ricardo/Data/GIAB/HG002_SV_benchmark/HG002_SVs_Tier1_v0.6.bed",
		ref = REFERENCE_FASTA,
		outdir = OUT_FOLDER + "/sv_genotyping/{sample}/bench/"
	conda:
		"../envs/truvari.yaml"
	shell:
		"rm -rf {params.outdir}; "
		#"zcat {input.vcf} | grep -E -v \"BREAKPOINT|COVERAGE\" | bcftools sort -Oz -o {input.vcf}.filt.vcf.gz; "
		"zcat {input.vcf} | grep -E \"^#|AGGREGATED\" | bcftools sort -Oz -o {input.vcf}.filt.vcf.gz; "
		"tabix -p vcf {input.vcf}.filt.vcf.gz; "
		"truvari bench \
 			-b {params.base} \
 			-c {input.vcf}.filt.vcf.gz \
 			-o {params.outdir} \
 			-f {params.ref} \
 			--includebed {params.bed} \
 			-p 0 \
 			-r 2000 \
 			--giabreport; "
