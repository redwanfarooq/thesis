##########################################################################################
# Snakemake rule for Cell Ranger ARC
# Author: Redwan Farooq
# Requires functions from resources/scripts/rule.py
# Requires outputs from resources/rules/bcl2fastq.smk
##########################################################################################

# Define rule
rule cellranger_arc:
	input: lambda wildcards: get_count_inputs(wildcards, lib_types={"GEX", "ATAC"}, info=info, read_trim=True if "read_trim" in config.keys() else False),
	output: os.path.abspath("stamps/cellranger_arc/{sample}.stamp")
	log: os.path.abspath("logs/cellranger_arc/{sample}.log")
	threads: 1
	params: 
		librarysheet_path = os.path.abspath(os.path.join(config.get("metadata_dir", "metadata"), "cellranger_arc")),
		reference = config["cellranger_reference"],
		custom_flags = config.get("cellranger_args", ""),
		output_path = os.path.join(config["output_dir"], "cellranger_arc")
	envmodules: "cellranger-arc/2.0.2"
	message: "Making GEX and ATAC count matrix for {wildcards.sample}"
	shell:
		"""
		( \
		mkdir -p stamps/cellranger_arc && \
		mkdir -p {params.output_path} && \
		cd {params.output_path} && \
		cellranger-arc count \
			--id={wildcards.sample} \
			--libraries={params.librarysheet_path}/{wildcards.sample}.csv \
			--reference={params.reference} \
			{params.custom_flags} \
			--localcores={threads} && \
		touch {output} \
		) > {log} 2>&1
		"""


# Set rule targets
cellranger_arc = [f"stamps/cellranger_arc/{sample}.stamp" for sample in samples]