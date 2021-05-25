rule delly_cnv:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai",
		bcf = OUT_FOLDER + "/sv_discovery/delly/{sample}/{sample}.bcf"
	output:
		vcf = OUT_FOLDER + "/sv_discovery/delly_cnv/{sample}/{sample}.delly_cnv.vcf.gz",
		tbi = OUT_FOLDER + "/sv_discovery/delly_cnv/{sample}/{sample}.delly_cnv.vcf.gz.tbi"
	conda:
		"../../envs/delly.yaml"
	params:
		delly_map = LIB_DIR + "/delly_maps/Homo_sapiens.GRCh" + REF_BUILD + ".dna.primary_assembly.fa.r101.s501.blacklist.gz",
		bcf = OUT_FOLDER + "/sv_discovery/delly_cnv/{sample}/{sample}.delly_cnv.bcf"
	shell:
		"delly cnv -g {REFERENCE_FASTA} -o {params.bcf} -m {params.delly_map} -l {input.bcf} $(readlink {input.bam}); "
		"bcftools sort -Oz -o {output.vcf} {params.bcf}; "
		"tabix -p vcf {output.vcf}; "