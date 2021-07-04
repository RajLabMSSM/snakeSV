rule smoove:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		outdir = temp(directory(OUT_FOLDER + "/sv_discovery/smoove/{sample}/outdir/")),
		vcf = OUT_FOLDER + "/sv_discovery/smoove/{sample}/{sample}.smoove.vcf.gz",
		tbi = OUT_FOLDER + "/sv_discovery/smoove/{sample}/{sample}.smoove.vcf.gz.tbi"
	params:
		exclude = LIB_DIR + "/smoove_regions/GRCh" + REF_BUILD + ".exclude.bed"
	conda:
		SNAKEDIR + "envs/smoove.yaml"
	shell:
		"smoove call --outdir {output.outdir} \
			--exclude {params.exclude} \
			--name {wildcards.sample} \
			--fasta {REFERENCE_FASTA} --genotype {input.bam}; "
		"zgrep -v 'SVTYPE=BND' {output.outdir}/{wildcards.sample}-smoove.genotyped.vcf.gz | \
			bcftools sort -Oz -o {output.vcf}; "
		"tabix -p vcf {output.vcf}; "
