###################################################################################################
## Other
###################################################################################################
rule somalier_extract:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		OUT_FOLDER + "/somalier_extracted/{sample}.somalier"
	params:
		ref = REF_BUILD
	shell:
		"{LIB_DIR}/somalier extract -d {OUT_FOLDER}/somalier_extracted/ \
			--sites {LIB_DIR}/somalier_sites/sites.GRCh{params.ref}.vcf.gz \
			-f {REFERENCE_FASTA} {input.bam}; "

rule somalier:
	input:
		expand(OUT_FOLDER + "/somalier_extracted/{sample}.somalier", sample=participant_id)
	output:
		OUT_FOLDER + "/somalier/somalier.somalier-ancestry.html"
	params:
		ref = REF_BUILD
	container:
		"docker://brentp/somalier"
	shell:
		"somalier relate -o {OUT_FOLDER}/somalier/somalier {input}; "
		"somalier ancestry --n-pcs=20 -o {OUT_FOLDER}/somalier/ \
			--labels {LIB_DIR}/somalier_sites/ancestry-labels-1kg.tsv \
			{LIB_DIR}/somalier_sites/1kg-somalier/*.somalier ++ {input}; "

rule svtk_collect_pers:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		splitfile = OUT_FOLDER + "/support/{sample}.splitfile",
		discfile = OUT_FOLDER + "/support/{sample}.discfile"
	conda:
		SNAKEDIR + "envs/svtk.yaml"	
	shell:
		"svtk collect-pesr {input.bam} {wildcards.sample} {output.splitfile} {output.discfile}; "
