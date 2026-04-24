##########################################################################################
# Snakemake rule for Cell Ranger Multi
# Author: Redwan Farooq
# Requires functions from resources/scripts/rule.py
# Requires outputs from resources/rules/bcl2fastq.smk
##########################################################################################

# Define rule
rule cellranger_multi:
	input: lambda wildcards: get_count_inputs(wildcards, lib_types={"GEX", "ADT", "HTO", "CRISPR", "BCR", "TCR"}, info=info, read_trim=True if "read_trim" in config.keys() else False),
	output: os.path.abspath("stamps/cellranger_multi/{sample}.stamp")
	log: os.path.abspath("logs/cellranger_multi/{sample}.log")
	threads: 1
	params: 
		librarysheet_path = os.path.abspath(os.path.join(config.get("metadata_dir", "metadata"), "cellranger")),
		custom_flags = config.get("cellranger_args", ""),
		output_path = os.path.join(config["output_dir"], "cellranger_multi")
	envmodules: "cellranger/9.0.0"
	message: "Running Cell Ranger Multi for {wildcards.sample}"
	shell:
		"""
		( \
		mkdir -p stamps/cellranger_multi && \
		mkdir -p {params.output_path} && \
		cd {params.output_path} && \
		cellranger multi \
			--id={wildcards.sample} \
			--csv={params.librarysheet_path}/{wildcards.sample}.csv \
			{params.custom_flags} \
			--localcores={threads} && \
		touch {output} \
		) > {log} 2>&1
		"""


# Set rule targets
cellranger_multi = [f"stamps/cellranger_multi/{sample}.stamp" for sample in samples]