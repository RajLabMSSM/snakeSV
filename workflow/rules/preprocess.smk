###################################################################################################
## Preprocess
###################################################################################################
rule symlink_files:
	input:
		config["SAMPLE_KEY"]
	output:
		expand(OUT_FOLDER + "/input/{sample}.bam", sample=participant_id)
	run:
		for label,row in f_input_list.iterrows():
			if os.path.isfile(OUT_FOLDER + "/input/" + row['participant_id'] + ".bam"):
				print("Link to BAM already created: " + row['participant_id'])
			else:
				os.system("ln -f -s $(realpath " + row['bam'] + ") " + OUT_FOLDER + "/input/" + row['participant_id'] + ".bam")

rule indexBAM:
	input:
		OUT_FOLDER + "/input/{sample}.bam"
	output:
		OUT_FOLDER + "/input/{sample}.bam.bai"
	conda:
		srcdir("../envs/bamtools.yaml")
	shell:
		"bamtools index -in {input}"
