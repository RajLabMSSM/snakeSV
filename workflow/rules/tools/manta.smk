rule manta:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		outdir = directory(OUT_FOLDER + "/sv_discovery/manta/{sample}/outdir/"),
		vcf = temp(OUT_FOLDER + "/sv_discovery/manta/{sample}/{sample}.manta.vcf")
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
		"zgrep -v 'SVTYPE=BND' {output.outdir}/results/variants/diploidSV.vcf.gz > {output.vcf}; "

rule compress_manta:
	input:
		vcf = OUT_FOLDER + "/sv_discovery/manta/{sample}/{sample}.manta.vcf"
	output:
		vcf = OUT_FOLDER + "/sv_discovery/manta/{sample}/{sample}.manta.vcf.gz",
		tbi = OUT_FOLDER + "/sv_discovery/manta/{sample}/{sample}.manta.vcf.gz.tbi"
	conda:
		SNAKEDIR + "envs/bcftools.yaml"
	shell:
		"bcftools sort -Oz -o {output.vcf} {input.vcf}; "
		"tabix -p vcf {output.vcf}; "

