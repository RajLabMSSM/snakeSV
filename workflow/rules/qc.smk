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
	log: 
		OUT_FOLDER + "/log/qc/{sample}/CollectMultipleMetrics.log"
	conda:
		SNAKEDIR + "envs/picard.yaml"
	shell:
		"picard CollectMultipleMetrics "
		"TMP_DIR={TMP_DIR} I={input.bam} AS=true R={REFERENCE_FASTA} VALIDATION_STRINGENCY=SILENT PROGRAM=CollectAlignmentSummaryMetrics "
		"PROGRAM=CollectInsertSizeMetrics OUTPUT={params} > {log} 2>&1"

rule flagstat:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		OUT_FOLDER + "/qc/{sample}/{sample}.flagstat"
	log: 
		OUT_FOLDER + "/log/qc/{sample}/flagstat.log"
	conda:
		SNAKEDIR + "envs/sambamba.yaml"
	shell:
		"sambamba view -h -f bam -F 'not secondary_alignment' {input.bam} | samtools flagstat /dev/stdin 1> {output}  2>{log}"

rule stats:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		OUT_FOLDER + "/qc/{sample}/{sample}.stats"
	log: 
		OUT_FOLDER + "/log/qc/{sample}/stats.log"
	conda:
		SNAKEDIR + "envs/sambamba.yaml"
	shell:
		"sambamba view -h -f bam -F 'not secondary_alignment' {input.bam} | bamtools stats -in /dev/stdin -insert 1> {output} 2>{log}"

rule wgs:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		OUT_FOLDER + "/qc/{sample}/{sample}.wgs"
	log: 
		OUT_FOLDER + "/log/qc/{sample}/wgs.log"		
	conda:
		SNAKEDIR + "envs/picard.yaml"
	shell:
		"picard CollectWgsMetrics TMP_DIR={TMP_DIR} I={input.bam} O={output} R={REFERENCE_FASTA} VALIDATION_STRINGENCY=SILENT > {log} 2>&1"

rule complexity:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		OUT_FOLDER + "/qc/{sample}/{sample}.complexity"
	log: 
		OUT_FOLDER + "/log/qc/{sample}/complexity.log"			
	conda:
		SNAKEDIR + "envs/picard.yaml"
	shell:
		"picard -Xms8g -Xmx16g EstimateLibraryComplexity TMP_DIR={TMP_DIR} I={input.bam} O={output} VALIDATION_STRINGENCY=SILENT > {log} 2>&1"

rule sexCheck:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		OUT_FOLDER + "/qc/{sample}/{sample}.sexCheck"
	params:
		script = SNAKEDIR + "/scripts/sexCheck.sh",
		DICT = config["DICT"]
	log: 
		OUT_FOLDER + "/log/qc/{sample}/sexCheck.log"		
	conda:
		SNAKEDIR + "envs/sambamba.yaml"
	shell:
		"sh {params.script} {wildcards.sample} {input.bam} {params.DICT} {OUT_FOLDER}/qc/{wildcards.sample} > {log} 2>&1"

rule aneuploidyCheck:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		OUT_FOLDER + "/qc/{sample}/{sample}.aneuploidyCheck"
	params:
		script = SNAKEDIR + "/scripts/chrCopyCount.sh",
		DICT = config["DICT"],
		NMASK = config["NMASK"]
	log: 
		OUT_FOLDER + "/log/qc/{sample}/aneuploidyCheck.log"			
	conda:
		SNAKEDIR + "envs/sambamba.yaml"
	shell:
		"sh {params.script} {wildcards.sample} {input.bam} {params.DICT} {OUT_FOLDER}/qc/{wildcards.sample} {params.NMASK} > {log} 2>&1"

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
		SNAKEDIR + "envs/multiqc.yaml"
	shell:
		"multiqc -f -o {OUT_FOLDER}/qc {OUT_FOLDER}/qc; "
