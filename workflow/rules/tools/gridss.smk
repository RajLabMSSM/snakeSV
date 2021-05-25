rule manta:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
        bai = OUT_FOLDER + "/input/{sample}.bam.bai"
	output:
		vcf = OUT_FOLDER + "/sv_discovery/manta/{sample}/{sample}.manta.vcf.gz",
		tbi = OUT_FOLDER + "/sv_discovery/manta/{sample}/{sample}.manta.vcf.gz.tbi"
    conda:
        "../../envs/gridss.yaml"
	params:
		ref = REF_BUILD
	shell:
		"conda deactivate; conda activate SV_pipeline_py27; "
		"configManta.py \
			--bam {input.bam} \
			--REFERENCE_FASTA {REFERENCE_FASTA} \
			--runDir {OUT_FOLDER}/sv_discovery/manta/{wildcards.sample}/; "
		"{OUT_FOLDER}/sv_discovery/manta/{wildcards.sample}/runWorkflow.py; "
		"bcftools sort -Oz -o {output.vcf} {OUT_FOLDER}/sv_discovery/manta/{wildcards.sample}/results/variants/diploidSV.vcf.gz; "
		"tabix -p vcf {output.vcf}; "

rule gridss_s:  # single-sample analysis
    input:
        fasta = get_fasta(),
        fai = get_faidx(),  # bwa index files also required
        bam = get_bam("{path}/{sample}"),
        bai = get_bai("{path}/{sample}")
    params:
        excl_opt = "BLACKLIST=" + get_bed() if exclude_regions() else ""
    output:
        os.path.join("{path}/{sample}", get_outdir("gridss"), "gridss" +
                     get_filext("vcf"))
    conda:
        "../environment.yaml"
    threads:
        get_nthreads("gridss")
    resources:
        mem_mb = get_memory("gridss"),
        tmp_mb = get_tmpspace("gridss")
    shell:
    	"gridss.sh --reference <reference.fa> \
    		--output <output.vcf.gz> \
    		--assembly <assembly.bam> \
    		[--threads n] [--jar gridss.jar] \
    		[--workingdir <directory>] \
    		[--jvmheap 30g] \
    		[--blacklist <blacklist.bed>] \
    		[--steps All|PreProcess|Assemble|Call] \
    		[--configuration gridss.properties] \
    		[--maxcoverage 50000] \
    		[--labels input1,input2,...] \
    		input1.bam [input2.bam [...]]"
        """

        set -x
        # if 'tmpspace' set to >0MB use TMPDIR otherwise use OUTDIR
        OUTDIR="$(dirname "{output}")"
        PREFIX="$(basename "{output}" .vcf)"
        OUTFILE="${{OUTDIR}}/${{PREFIX}}.unfiltered.vcf"
        TMP=$([ "{resources.tmp_mb}" -eq "0" ] &&
            echo "${{OUTDIR}}" || echo "${{TMPDIR}}")
        # set JVM max. heap size dynamically (in GB)
        # N.B. don't allocate >31G due to Compressed Oops and JDK-8029679
        MAX_HEAP=$(LC_ALL=C printf "%.f" $(bc <<< "scale=2; \
            {resources.mem_mb} / 1024 * .8")) # max. 80% of requested mem
        MAX_HEAP=$([ "${{MAX_HEAP}}" -gt "31" ] && echo "31g" ||
            echo "${{MAX_HEAP}}g")
        export _JAVA_OPTIONS="-Djava.io.tmpdir=${{TMP}} -Xmx${{MAX_HEAP}}"
        # run dummy or real job
        if [ "{config[echo_run]}" -eq "1" ]; then
            echo "{input}" "${{TMP}}" > "{output}"
        else
            # clean-up outdir prior to SV calling
            rm -fr ${{OUTDIR}}/*gridss* &&
            gridss gridss.CallVariants \
                WORKER_THREADS={threads} \
                REFERENCE_SEQUENCE="{input.fasta}" \
                {params.excl_opt} \
                INPUT="{input.bam}" \
                OUTPUT="${{OUTFILE}}" \
                ASSEMBLY="${{OUTDIR}}/gridss_assembly.bam" \
                WORKING_DIR="${{TMP}}" \
                TMP_DIR="${{TMP}}/gridss.${{RANDOM}}" &&
            # SV quality filtering
            bcftools filter \
                -O v `# uncompressed VCF format` \
                -o "{output}" \
                -i "FILTER == '.'" \
                "${{OUTFILE}}"
        fi
        """