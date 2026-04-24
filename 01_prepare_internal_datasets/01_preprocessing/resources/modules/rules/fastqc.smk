##########################################################################################
# Snakemake rule for FastQC
# Author: Redwan Farooq
# Requires functions from resources/scripts/rule.py
# Requires outputs from resources/rules/bcl2fastq.smk or resources/rules/trimfastq.smk
##########################################################################################

# Define rule
rule fastqc:
	input: lambda wildcards: get_fastqc_inputs(wildcards, info=info, read_trim=True if "read_trim" in config.keys() else False)
	output: os.path.abspath("stamps/fastqc/{lib}.stamp")
	log: os.path.abspath("logs/fastqc/{lib}.log")
	threads: 1
	params:
		fastqs = lambda wildcards: get_fastqc_fastqs(wildcards, info=info, output_dir=config["output_dir"]),
		custom_flags = config.get("fastqc_args", ""),
		output_path = os.path.join(config["output_dir"], "qc/fastqc")
	# conda: "fastqc"
	envmodules: "fastqc/0.11.9"
	message: "Performing QC for FASTQ files for {wildcards.lib}"
	shell:
		"""
		( \
		mkdir -p stamps/fastqc && \
		mkdir -p {params.output_path}/{wildcards.lib} && \
		fastqc \
			--threads={threads} \
			{params.custom_flags} \
			--outdir={params.output_path}/{wildcards.lib} \
			{params.fastqs} && \
		touch {output} \
		) > {log} 2>&1
		"""

# Set rule targets
fastqc = [f"stamps/fastqc/{lib}.stamp" for lib in libs]