##########################################################################################
# Snakemake rule for MultiQC
# Author: Redwan Farooq
# Requires outputs from resources/rules/fastqc.smk
##########################################################################################

# Define rule
rule multiqc:
	input: [os.path.abspath(path) for path in expand("stamps/fastqc/{lib}.stamp", lib=libs.keys())]
	output: os.path.abspath("stamps/multiqc/multiqc.stamp")
	log: os.path.abspath("logs/multiqc/multiqc.log")
	threads: 1
	params:
		input_path = os.path.join(config["output_dir"], "qc"),
		custom_flags = config.get("multiqc_args", ""),
		output_path = os.path.join(config["output_dir"], "qc/multiqc") # DO NOT CHANGE - downstream rules will search for summary statistics in this directory
	conda: "multiqc"
	# envmodules: "multiqc/1.14"
	message: "Aggregating QC reports"
	shell:
		"""
		( \
		mkdir -p stamps/multiqc && \
		mkdir -p {params.output_path} && \
		multiqc \
			--outdir={params.output_path} \
			{params.custom_flags} \
			{params.input_path} && \
		touch {output} \
		) > {log} 2>&1
		"""


# Set rule targets
multiqc = ["stamps/multiqc/multiqc.stamp"]