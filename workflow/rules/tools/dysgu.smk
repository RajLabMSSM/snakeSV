rule dysgu:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		tmp_dir = temp(directory(OUT_FOLDER + "/sv_discovery/dysgu/{sample}/tmp")),
		vcf = temp(OUT_FOLDER + "/sv_discovery/dysgu/{sample}/{sample}.dysgu.vcf")
	conda:
		SNAKEDIR + "envs/dysgu.yaml"
	shell:
		"dysgu run {REFERENCE_FASTA} {output.tmp_dir} {input.bam} > {output.vcf}; "

rule compress_dysgu:
	input:
		vcf = OUT_FOLDER + "/sv_discovery/dysgu/{sample}/{sample}.dysgu.vcf"
	output:
		vcf = OUT_FOLDER + "/sv_discovery/dysgu/{sample}/{sample}.dysgu.vcf.gz",
		tbi = OUT_FOLDER + "/sv_discovery/dysgu/{sample}/{sample}.dysgu.vcf.gz.tbi"
	conda:
		SNAKEDIR + "envs/bcftools.yaml"
	shell:
		"bcftools sort -Oz -o {output.vcf} {input.vcf}; "
		"tabix -p vcf {output.vcf}; "
		
