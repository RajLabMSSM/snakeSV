###################################################################################################
## Other
###################################################################################################
if( "ANNOTATION_BED" not in config ):
	rule annotate_bed:
		input:
			vcf = OUT_FOLDER + "/merged_cohort/gt_merged.vcf",
			vcfgz = OUT_FOLDER + "/merged_cohort/gt_merged.vcf.gz",
			tbi = OUT_FOLDER + "/merged_cohort/gt_merged.vcf.gz.tbi"
		output:
			vcf = OUT_FOLDER + "/merged_cohort/gt_merged.annot.vcf.gz",
			tbi = OUT_FOLDER + "/merged_cohort/gt_merged.annot.vcf.gz.tbi"
		conda:
			"../envs/svtk.yaml"		
		params:
			gencode = config["GENCODE_GTF"]
		shell:
			"svtk annotate \
				--gencode {params.gencode} \
				{input.vcfgz} {output.vcf}; "
			"tabix -p vcf {output.vcf}; "
else:
	rule annotate_bed:
		input:
			vcf = OUT_FOLDER + "/merged_cohort/gt_merged.vcf",
			vcfgz = OUT_FOLDER + "/merged_cohort/gt_merged.vcf.gz",
			tbi = OUT_FOLDER + "/merged_cohort/gt_merged.vcf.gz.tbi"
		output:
			vcf = OUT_FOLDER + "/merged_cohort/gt_merged.annot.vcf.gz",
			tbi = OUT_FOLDER + "/merged_cohort/gt_merged.annot.vcf.gz.tbi"
		conda:
			"../envs/svtk.yaml"		
		params:
			bed = config["ANNOTATION_BED"],
			bed_cat = OUT_FOLDER + "/merged_cohort/annot.bed",
			gencode = config["GENCODE_GTF"]
		shell:
			"cat {params.bed} > {params.bed_cat}; "
			"svtk annotate \
				--gencode {params.gencode} \
				--noncoding {params.bed_cat} \
				{input.vcfgz} {output.vcf}; "
			"tabix -p vcf {output.vcf}; "
