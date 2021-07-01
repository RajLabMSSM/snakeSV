rule delly:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		vcf = OUT_FOLDER + "/sv_discovery/delly/{sample}/{sample}.delly.vcf.gz",
		tbi = OUT_FOLDER + "/sv_discovery/delly/{sample}/{sample}.delly.vcf.gz.tbi",
		bcf = OUT_FOLDER + "/sv_discovery/delly/{sample}/{sample}.bcf"
	conda:
		SNAKEDIR + "envs/delly.yaml"
	params:
		exclude = LIB_DIR + "/delly_maps/human.hg" + REF_BUILD + ".excl.tsv",
	shell:
		"delly call -g {REFERENCE_FASTA} -o {output.bcf} -x {params.exclude} $(readlink {input.bam}); "
		"bcftools sort -Oz -o {output.vcf} {output.bcf}; "
		"tabix -p vcf {output.vcf}; "
