rule manta:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		outdir = directory(OUT_FOLDER + "/sv_discovery/manta/{sample}/outdir/"),
		vcf = OUT_FOLDER + "/sv_discovery/manta/{sample}/{sample}.manta.vcf.gz",
		tbi = OUT_FOLDER + "/sv_discovery/manta/{sample}/{sample}.manta.vcf.gz.tbi"
	conda:
		SNAKEDIR + "envs/manta.yaml"
	params:
		ref = REF_BUILD
	shell:
		"rm -rf {OUT_FOLDER}/sv_discovery/manta/{wildcards.sample}/; "
		"configManta.py \
			--bam {input.bam} \
			--referenceFasta {REFERENCE_FASTA} \
			--runDir {output.outdir}; "
		"{output.outdir}/runWorkflow.py; "
		"zgrep -v 'SVTYPE=BND' {output.outdir}/results/variants/diploidSV.vcf.gz | \
			bcftools sort -Oz -o {output.vcf}; "
		"tabix -p vcf {output.vcf}; "
