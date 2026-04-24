##########################################################################################
# Snakemake rule for macs2
# Author: Redwan Farooq
# Requires outputs from resources/rules/chromap.smk
##########################################################################################

scripts_dir = config.get("scripts_dir", "resources/scripts")

# Define rule
rule macs2:
	input: os.path.abspath("stamps/chromap/{sample}.stamp")
	output: os.path.abspath("stamps/macs2/{sample}.stamp")
	log: os.path.abspath("logs/macs2/{sample}.log")
	threads: 1
	params:
		custom_flags = config.get("macs2_args", ""),
		script_path = scripts_dir if os.path.isabs(scripts_dir) else os.path.join(workflow.basedir, scripts_dir),
		output_path = os.path.join(config["output_dir"], "chromap_macs2") # DO NOT CHANGE - downstream rules will search for mapping statistics in this directory
	conda: "macs2"
	# envmodules:
	# 	"MACS2/2.2.9.1",
	# 	"R-cbrg"
	message: "Making ATAC count matrix for {wildcards.sample}"
	shell:
		"""
		( \
		mkdir -p stamps/macs2 && \
		mkdir -p {params.output_path}/{wildcards.sample} && \
		macs2 callpeak \
			-t {params.output_path}/{wildcards.sample}/fragments.tsv.gz \
			-n macs2 \
			--outdir {params.output_path}/{wildcards.sample} \
			{params.custom_flags} && \
		{params.script_path}/count_peaks.R \
			--fragments {params.output_path}/{wildcards.sample}/fragments.tsv.gz \
			--peaks {params.output_path}/{wildcards.sample}/macs2_peaks.narrowPeak \
			--threads {threads} && \
		find {params.output_path}/{wildcards.sample}/raw_feature_bc_matrix -type f -name 'barcodes.tsv' -exec gzip -f {{}} \; && \
		find {params.output_path}/{wildcards.sample}/raw_feature_bc_matrix -type f -name 'features.tsv' -exec gzip -f {{}} \; && \
		find {params.output_path}/{wildcards.sample}/raw_feature_bc_matrix -type f -name 'matrix.mtx' -exec gzip -f {{}} \; && \
		touch {output} \
		) > {log} 2>&1
		"""


# Set rule targets
macs2 = [f"stamps/macs2/{sample}.stamp" for sample in samples]