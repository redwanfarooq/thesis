##########################################################################################
# Snakemake rule for Cell Ranger ATAC
# Author: Redwan Farooq
# Requires functions from resources/scripts/rule.py
# Requires outputs from resources/rules/bcl2fastq.smk
##########################################################################################

# Define rule
rule cellranger_atac:
	input: lambda wildcards: get_count_inputs(wildcards, lib_types={"ATAC"}, info=info, read_trim=True if "read_trim" in config.keys() else False),
	output: os.path.abspath("stamps/cellranger_atac/{sample}.stamp")
	log: os.path.abspath("logs/cellranger_atac/{sample}.log")
	threads: 1
	params: 
		fastqdirs = lambda wildcards: get_count_fastqdirs(wildcards, lib_types={"ATAC"}, info=info, output_dir=config["output_dir"]),
		reference = config["cellranger_reference"],
		custom_flags = config.get("cellranger_args", ""),
		output_path = os.path.join(config["output_dir"], "cellranger_atac")
	envmodules: "cellranger-atac/2.1.0"
	message: "Making ATAC count matrix for {wildcards.sample}"
	shell:
		"""
		( \
		mkdir -p stamps/cellranger_atac && \
		mkdir -p {params.output_path} && \
		cd {params.output_path} && \
		cellranger-atac count \
			--id={wildcards.sample} \
			--sample={wildcards.sample} \
			--fastqs={params.fastqdirs} \
			--reference={params.reference} \
			{params.custom_flags} \
			--localcores={threads} && \
		touch {output} \
		) > {log} 2>&1
		"""


# Set rule targets
cellranger_atac = [f"stamps/cellranger_atac/{sample}.stamp" for sample in samples]