rule smoove:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		outdir = temp(directory(OUT_FOLDER + "/sv_discovery/smoove/{sample}/outdir/")),
		vcf = temp(OUT_FOLDER + "/sv_discovery/smoove/{sample}/{sample}.smoove.vcf")
	conda:
		SNAKEDIR + "envs/smoove.yaml"
	shell:
		"smoove call --outdir {output.outdir} \
			--name {wildcards.sample} \
			--fasta {REFERENCE_FASTA} --genotype {input.bam}; "
		"zgrep -v 'SVTYPE=BND' {output.outdir}/{wildcards.sample}-smoove.genotyped.vcf.gz > {output.vcf}; "

rule compress_smoove:
	input:
		vcf = OUT_FOLDER + "/sv_discovery/smoove/{sample}/{sample}.smoove.vcf"
	output:
		vcf = OUT_FOLDER + "/sv_discovery/smoove/{sample}/{sample}.smoove.vcf.gz",
		tbi = OUT_FOLDER + "/sv_discovery/smoove/{sample}/{sample}.smoove.vcf.gz.tbi"
	conda:
		SNAKEDIR + "envs/bcftools.yaml"
	shell:
		"bcftools sort -Oz -o {output.vcf} {input.vcf}; "
		"tabix -p vcf {output.vcf}; "
