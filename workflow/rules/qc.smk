###################################################################################################
## Sample QC
###################################################################################################
rule CollectMultipleMetrics:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		OUT_FOLDER + "/qc/{sample}/{sample}.metrics.insert_size_metrics"
	params:
		OUT_FOLDER + "/qc/{sample}/{sample}.metrics"
	conda:
		srcdir("../envs/picard.yaml")
	shell:
		"picard CollectMultipleMetrics "
		"TMP_DIR={TMP_DIR} I={input.bam} AS=true R={REFERENCE_FASTA} VALIDATION_STRINGENCY=SILENT PROGRAM=CollectAlignmentSummaryMetrics "
		"PROGRAM=CollectInsertSizeMetrics OUTPUT={params}"

rule flagstat:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		OUT_FOLDER + "/qc/{sample}/{sample}.flagstat"
	conda:
		srcdir("../envs/sambamba.yaml")
	shell:
		"sambamba view -h -f bam -F 'not secondary_alignment' {input.bam} | samtools flagstat /dev/stdin > {output}"

rule stats:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		OUT_FOLDER + "/qc/{sample}/{sample}.stats"
	conda:
		srcdir("../envs/sambamba.yaml")
	shell:
		"sambamba view -h -f bam -F 'not secondary_alignment' {input.bam} | bamtools stats -in /dev/stdin -insert > {output}"

rule wgs:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		OUT_FOLDER + "/qc/{sample}/{sample}.wgs"
	conda:
		srcdir("../envs/picard.yaml")
	shell:
		"picard CollectWgsMetrics TMP_DIR={TMP_DIR} I={input.bam} O={output} R={REFERENCE_FASTA} VALIDATION_STRINGENCY=SILENT"

rule complexity:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		OUT_FOLDER + "/qc/{sample}/{sample}.complexity"
	conda:
		srcdir("../envs/picard.yaml")
	shell:
		"picard EstimateLibraryComplexity TMP_DIR={TMP_DIR} I={input.bam} O={output} VALIDATION_STRINGENCY=SILENT"

rule sexCheck:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		OUT_FOLDER + "/qc/{sample}/{sample}.sexCheck"
	conda:
		srcdir("../envs/sambamba.yaml")
	params:
		script = "workflow/scripts/sexCheck.sh",
		DICT = config["DICT"]
	shell:
		"sh {params.script} {wildcards.sample} {input.bam} {params.DICT} {OUT_FOLDER}/qc/{wildcards.sample}"

rule aneuploidyCheck:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		OUT_FOLDER + "/qc/{sample}/{sample}.aneuploidyCheck"
	conda:
		srcdir("../envs/sambamba.yaml")
	params:
		script = "workflow/scripts/chrCopyCount.sh",
		DICT = config["DICT"],
		NMASK = config["NMASK"]
	shell:
		"sh {params.script} {wildcards.sample} {input.bam} {params.DICT} {OUT_FOLDER}/qc/{wildcards.sample} {params.NMASK}"

rule multiqc:
	input:
		expand(OUT_FOLDER + "/qc/{sample}/{sample}.metrics.insert_size_metrics", sample=participant_id),
		expand(OUT_FOLDER + "/qc/{sample}/{sample}.flagstat", sample=participant_id),
		expand(OUT_FOLDER + "/qc/{sample}/{sample}.stats", sample=participant_id),
		expand(OUT_FOLDER + "/qc/{sample}/{sample}.wgs", sample=participant_id),
		expand(OUT_FOLDER + "/qc/{sample}/{sample}.complexity", sample=participant_id),
		expand(OUT_FOLDER + "/qc/{sample}/{sample}.sexCheck", sample=participant_id),
		expand(OUT_FOLDER + "/qc/{sample}/{sample}.aneuploidyCheck", sample=participant_id)
	output:
		OUT_FOLDER + "/qc/multiqc_report.html"
	conda:
		srcdir("../envs/multiqc.yaml")
	shell:
		"multiqc -f -o {OUT_FOLDER}/qc {OUT_FOLDER}/qc; "
