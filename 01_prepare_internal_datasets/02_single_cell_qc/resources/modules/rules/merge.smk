##########################################################################################
# Snakemake rule for merging multimodal count matrices
# Author: Redwan Farooq
# Requires functions from resources/scripts/rule.py
##########################################################################################

scripts_dir = config.get("scripts_dir", "resources/scripts")

# Define rule
rule merge:
	output: os.path.join(config["output_dir"]["data"], f"merge/{{sample}}/{'filtered' if config.get('pre_filtered', None) and 'cellbender' not in module_rules else 'raw'}_feature_bc_matrix.h5")
	log: os.path.abspath("logs/merge/{sample}.log")
	threads: 1
	params:
		script_path = scripts_dir if os.path.isabs(scripts_dir) else os.path.join(workflow.basedir, scripts_dir),
		flags = lambda wildcards: get_merge_flags(wildcards, gex=config.get("gex_matrix", None), atac=config.get("atac_matrix", None), adt=config.get("adt_matrix", None), adt_prefix=config.get("adt_prefix", None))
	conda: "single_cell_qc"
	# envmodules: "R-cbrg"
	message: "Merging multimodal count matrices for {wildcards.sample}"
	shell:
		"""
		( \
		cd {params.script_path} && \
		./merge.R \
			{params.flags} \
			--output {output} \
		) > {log} 2>&1
		"""

	
# Set rule targets
merge = expand(os.path.join(config["output_dir"]["data"], f"merge/{{sample}}/{'filtered' if config.get('pre_filtered', None) and 'cellbender' not in module_rules else 'raw'}_feature_bc_matrix.h5"), sample=samples)