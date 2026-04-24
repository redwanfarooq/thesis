##########################################################################################
# Snakemake rule for filtering merged multimodal count matrices by cell barcode
# Author: Redwan Farooq
# Requires functions from resources/scripts/rule.py
# Requires outputs from resources/rules/libraries_qc.smk
##########################################################################################

scripts_dir = config.get("scripts_dir", "resources/scripts")

# Define rule
rule filter_barcodes:
	input: os.path.join(config["output_dir"]["reports"], "libraries_qc/{sample}.html")
	output: os.path.join(config["output_dir"]["data"], "finalised/{sample}/decontaminated_matrix.h5")
	log: os.path.abspath("logs/filter_barcodes/{sample}.log")
	threads: 1
	params:
		script_path = scripts_dir if os.path.isabs(scripts_dir) else os.path.join(workflow.basedir, scripts_dir),
		input_path = config["output_dir"]["data"],
		features_matrix = lambda wildcards: get_features_matrix(wildcards, data_dir=config["output_dir"]["data"], cellbender="cellbender" in module_rules, filtered=config.get("pre_filtered", None))
	conda: "single_cell_qc"
	# envmodules: "R-cbrg"
	message: "Filtering multimodal count matrices for {wildcards.sample}"
	shell:
		"""
		( \
		cd {params.script_path} && \
		./filter_barcodes.R \
			--input {params.features_matrix} \
			--barcodes {params.input_path}/libraries_qc/{wildcards.sample}/cell_barcodes.txt.gz \
			--output {output} \
		) > {log} 2>&1
		"""

	
# Set rule targets
filter_barcodes = expand(os.path.join(config["output_dir"]["data"], "finalised/{sample}/decontaminated_matrix.h5"), sample=samples)