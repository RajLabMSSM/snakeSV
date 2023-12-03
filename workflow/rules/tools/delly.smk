rule delly:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		bcf = OUT_FOLDER + "/sv_discovery/delly/{sample}/{sample}.bcf"
	conda:
		SNAKEDIR + "envs/delly.yaml"
	shell:
		"delly call -g {REFERENCE_FASTA} -o {output.bcf} $(readlink {input.bam}); "

rule compress_delly:
	input:
		bcf = OUT_FOLDER + "/sv_discovery/delly/{sample}/{sample}.bcf"
	output:
		vcf = OUT_FOLDER + "/sv_discovery/delly/{sample}/{sample}.delly.vcf.gz",
		tbi = OUT_FOLDER + "/sv_discovery/delly/{sample}/{sample}.delly.vcf.gz.tbi"
	conda:
		SNAKEDIR + "envs/bcftools.yaml"
	shell:
		"bcftools sort -Oz -o {output.vcf} {input.bcf}; "
		"tabix -p vcf {output.vcf}; "

