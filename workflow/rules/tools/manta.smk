rule manta:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		vcf = OUT_FOLDER + "/sv_discovery/manta/{sample}/{sample}.manta.vcf.gz",
		tbi = OUT_FOLDER + "/sv_discovery/manta/{sample}/{sample}.manta.vcf.gz.tbi"
	conda:
		"../../envs/manta.yaml"
	params:
		ref = REF_BUILD
	shell:
		"rm -rf {OUT_FOLDER}/sv_discovery/manta/{wildcards.sample}/; "
		"configManta.py \
			--bam {input.bam} \
			--referenceFasta {REFERENCE_FASTA} \
			--runDir {OUT_FOLDER}/sv_discovery/manta/{wildcards.sample}/; "
		"{OUT_FOLDER}/sv_discovery/manta/{wildcards.sample}/runWorkflow.py; "
		"zgrep -v 'SVTYPE=BND' {OUT_FOLDER}/sv_discovery/manta/{wildcards.sample}/results/variants/diploidSV.vcf.gz | \
			bcftools sort -Oz -o {output.vcf}; "
		"tabix -p vcf {output.vcf}; "
		