rule smoove:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		vcf = OUT_FOLDER + "/sv_discovery/smoove/{sample}/{sample}.smoove.vcf.gz",
		tbi = OUT_FOLDER + "/sv_discovery/smoove/{sample}/{sample}.smoove.vcf.gz.tbi"
	params:
		exclude = LIB_DIR + "/smoove_regions/GRCh" + REF_BUILD + ".exclude.bed"
	conda:
		SNAKEDIR + "envs/smoove.yaml"
	shell:
		"smoove call --outdir {OUT_FOLDER}/sv_discovery/smoove/ \
			--exclude {params.exclude} \
			--name {wildcards.sample} \
			--fasta {REFERENCE_FASTA} --genotype {input.bam}; "
		"zgrep -v 'SVTYPE=BND' {OUT_FOLDER}/sv_discovery/smoove/{wildcards.sample}-smoove.genotyped.vcf.gz | \
			bcftools sort -Oz -o {output.vcf}; "
		"tabix -p vcf {output.vcf}; "
