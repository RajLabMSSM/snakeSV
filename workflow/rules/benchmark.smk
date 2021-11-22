###################################################################################################
## Benchmarking against GiaB
# Deactivated in the pipeline.
###################################################################################################

# Benchmark raw SV callers
rule benchmark_raw:
	input:
		vcf = OUT_FOLDER + "/sv_discovery/{tool}/{sample}/{sample}.{tool}.vcf.gz"
	output:
		OUT_FOLDER + "/bench/{sample}/sv_discovery/{tool}/giab_report.txt"
	conda:
		SNAKEDIR + "envs/truvari.yaml"		
	params:
		base = GIAB_VCF,
		bed = GIAB_BED,
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

# Benchmark each SV caller after genotyping (e.g. Manta + GraphTyper, or Delly + GraphTyper)
rule benchmark_raw_gt:
	input:
		vcf = OUT_FOLDER + "/sv_discovery/{tool}/{sample}/{sample}.{tool}.gt.vcf.gz"
	output:
		OUT_FOLDER + "/bench/{sample}/sv_discovery/{tool}.gt/giab_report.txt"
	params:
		base = GIAB_VCF,
		bed = GIAB_BED,
		ref = REFERENCE_FASTA,
		outdir = OUT_FOLDER + "/bench/{sample}/sv_discovery/{tool}.gt/"
	conda:
		SNAKEDIR + "envs/truvari.yaml"
	shell:
		"rm -rf {params.outdir}; "
		"bcftools view -e QUAL==0 -Oz -o {input.vcf}.filt.vcf.gz {input.vcf}; "
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

# Benchmark each SV caller after joint-genotying. This may also include SV reference panels, if specified.
rule benchmark_gt:
	input:
		vcf = OUT_FOLDER + "/sv_genotyping/{sample}/{sample}.vcf.gz",
		filt_vcfgz = OUT_FOLDER + "/sv_genotyping/{sample}/{sample}_filt.vcf.gz"
	output:
		OUT_FOLDER + "/bench/{sample}/sv_genotyping/giab_report.txt"
	params:
		base = GIAB_VCF,
		bed = GIAB_BED,
		ref = REFERENCE_FASTA,
		outdir = OUT_FOLDER + "/bench/{sample}/sv_genotyping/"
	conda:
		SNAKEDIR + "envs/truvari.yaml"
	shell:
		"rm -rf {params.outdir}; "
		#"bcftools view -e QUAL==0 {input.vcf} | bcftools sort -Oz -o {input.vcf}.filt.vcf.gz; "
		#"tabix -p vcf {input.vcf}.filt.vcf.gz; "
		"truvari bench \
 			-b {params.base} \
 			-c {input.filt_vcfgz} \
 			-o {params.outdir} \
 			-f {params.ref} \
 			--includebed {params.bed} \
 			-p 0 \
 			-r 2000 \
 			--giabreport; "

# Benchmark SVs from the final genotyped and merged call set (collect only SVs present in the respective sample)
rule benchmark_merged_gt:
	input:
		vcf = OUT_FOLDER + "/merged_cohort/gt_merged.vcf.gz"
	output:
		vcf = OUT_FOLDER + "/bench/{sample}/merged/vcf.gz",
		report = OUT_FOLDER + "/bench/{sample}/merged/truvari/giab_report.txt"
	params:
		base = GIAB_VCF,
		bed = GIAB_BED,
		ref = REFERENCE_FASTA,
		outdir = OUT_FOLDER + "/bench/{sample}/merged/"
	conda:
		SNAKEDIR + "envs/truvari.yaml"
	shell:
		"mkdir -p {params.outdir}; "
		"bcftools view -c1 -s {wildcards.sample} {input.vcf} | bcftools sort -Oz -o {output.vcf}; "
		"tabix -p vcf {output.vcf}; "
		"rm -rf {params.outdir}/truvari; "
		"truvari bench \
 			-b {params.base} \
 			-c {output.vcf} \
 			-o {params.outdir}/truvari \
 			-f {params.ref} \
 			--includebed {params.bed} \
 			-p 0 \
 			-r 2000 \
 			--giabreport; "


