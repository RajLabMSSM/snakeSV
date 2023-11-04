###################################################################################################
## Annotations
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
			SNAKEDIR + "envs/svtk.yaml"		
		params:
			gencode = config["GENCODE_GTF"]
		priority: 0
		shell:
			"svtk annotate \
				--gencode {params.gencode} \
				{input.vcfgz} {output.vcf}; "
			"tabix -p vcf {output.vcf}; "
else: 
	rule annotate_bed:
		input:
			vcf = OUT_FOLDER + "/merged_cohort/gt_merged_full.vcf",
			vcfgz = OUT_FOLDER + "/merged_cohort/gt_merged_full.vcf.gz",
			tbi = OUT_FOLDER + "/merged_cohort/gt_merged_full.vcf.gz.tbi"
		output:
			tmp = temp(directory(OUT_FOLDER + "/merged_cohort/tmp")),
			vcf = OUT_FOLDER + "/merged_cohort/gt_merged.annot.vcf.gz",
			tbi = OUT_FOLDER + "/merged_cohort/gt_merged.annot.vcf.gz.tbi"
		conda:
			SNAKEDIR + "envs/svtk.yaml"
		params:
			bed = config["ANNOTATION_BED"],
			bed_cat = OUT_FOLDER + "/merged_cohort/annot.bed",
			gencode = config["GENCODE_GTF"],
			gencode_copy = OUT_FOLDER + "/merged_cohort/tmp/gencode.gtf",
			header = OUT_FOLDER + "/merged_cohort/tmp/header",
			variants = OUT_FOLDER + "/merged_cohort/tmp/variants"
		priority: 0
		shell:
			# svtk seems to strugle to process large files. 
			# breaking in small chunks for annotation. 
			"cat {params.bed} > {params.bed_cat}; "
			"mkdir -p {output.tmp}; "
			"cp {params.gencode} {params.gencode_copy}; "
			#grab the header
			"head -n 10000 {input.vcf} | grep '^#' > {params.header}; "
			#grab the non header lines
			"grep -v '^#' {input.vcf} > {params.variants}; "
			"curDir=$PWD; "
			#split into chunks with 100 lines
			"cd {output.tmp}; "
			"split -l 3000 variants; "
			#reattach the header to each and clean up
			"for i in x*; do "
				"cat header $i > $i.vcf && rm -f $i;"
				"svtk annotate \
					--gencode gencode.gtf \
					--noncoding ../annot.bed \
					$i.vcf $i.annot.vcf; "
			"done; "
			"cd $curDir; "
			"bcftools concat {output.tmp}/*.annot.vcf -Oz -o {output.vcf}; "
			"tabix -p vcf {output.vcf}; "
