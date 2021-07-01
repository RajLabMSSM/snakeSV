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
		SNAKEDIR + "envs/bamtools.yaml"
	shell:
		"bamtools index -in {input}"

rule indexREF:
	input:
		REFERENCE_FASTA
	output:
		REFERENCE_FASTA + ".fai"
	conda:
		SNAKEDIR + "envs/samtools.yaml"
	shell:
		"samtools faidx {input}"

rule makeDICT:
	input:
		REFERENCE_FASTA
	output:
		OUT_FOLDER + "/ref/" + config["REF_BUILD"] + ".dict"
	log: 
		OUT_FOLDER + "/log/qc/CreateSequenceDictionary.log"
	conda:
		SNAKEDIR + "envs/picard.yaml"
	shell:
		"picard CreateSequenceDictionary R={input} O={output} > {log} 2>&1"

