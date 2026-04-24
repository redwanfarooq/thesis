##########################################################################################
# Snakemake rule for cellbender
# Author: Redwan Farooq
# Requires functions from resources/scripts/rule.py
# Requires outputs from resources/rules/merge.smk
##########################################################################################

# Define rule
rule cellbender:
	input: os.path.join(config["output_dir"]["data"], "merge/{sample}/raw_feature_bc_matrix.h5")
	output: os.path.join(config["output_dir"]["reports"], "cellbender/{sample}.html")
	log: os.path.abspath("logs/cellbender/{sample}.log")
	threads: 1
	params:
		expected_cells_flag = lambda wildcards: get_expected_cells_flag(wildcards, info=info),
		ignore_features_flag = lambda wildcards, input: get_ignore_features_flag(hdf5=str(input[0]), regex=config.get("adt_isotype", None)),
		custom_flags = config.get("cellbender_args", ""),
		output_path = lambda wildcards: os.path.join(config["output_dir"]["data"], "cellbender") # DO NOT CHANGE - downstream rules will search for HDF5 files in this directory
	conda: "cellbender"
	# envmodules: "py-cellbender/0.3.0"
	message: "Running CellBender for {wildcards.sample}"
	shell:
		"""
		( \
		mkdir -p {params.output_path}/{wildcards.sample} && \
		cd {params.output_path}/{wildcards.sample} && \
		cellbender remove-background \
			--input {input} \
			--cpu-threads {threads} \
			{params.expected_cells_flag} \
			{params.ignore_features_flag} \
			{params.custom_flags} \
			--output=cellbender.h5 && \
		mv cellbender.h5 raw_feature_bc_matrix.h5 && \
		mv cellbender_filtered.h5 filtered_feature_bc_matrix.h5 && \
		mv cellbender_posterior.h5 posterior.h5 && \
		mv cellbender_cell_barcodes.csv cell_barcodes.csv && \
		mv cellbender_report.html {output} \
		) > {log} 2>&1
		"""

	
# Set rule targets
cellbender = expand(os.path.join(config["output_dir"]["reports"], "cellbender/{sample}.html"), sample=samples)