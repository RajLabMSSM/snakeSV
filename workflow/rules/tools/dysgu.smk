rule dysgu:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		tmp_dir = temp(directory(OUT_FOLDER + "/sv_discovery/dysgu/{sample}/tmp")),
		vcf = OUT_FOLDER + "/sv_discovery/dysgu/{sample}/{sample}.dysgu.vcf.gz",
		tbi = OUT_FOLDER + "/sv_discovery/dysgu/{sample}/{sample}.dysgu.vcf.gz.tbi"
	conda:
		SNAKEDIR + "envs/dysgu.yaml"
	shell:
		"dysgu run {REFERENCE_FASTA} {output.tmp_dir} {input.bam} | bcftools sort -Oz -o {output.vcf}; "
		"tabix -p vcf {output.vcf}; "
		