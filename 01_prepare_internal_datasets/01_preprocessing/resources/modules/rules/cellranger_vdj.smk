##########################################################################################
# Snakemake rule for Cell Ranger VDJ
# Author: Redwan Farooq
# Requires functions from resources/scripts/rule.py
# Requires outputs from resources/rules/bcl2fastq.smk
##########################################################################################

# Define rule
rule cellranger_vdj:
	input: lambda wildcards: get_count_inputs(wildcards, lib_types={"BCR", "TCR"}, info=info, read_trim=True if "read_trim" in config.keys() else False),
	output: os.path.abspath("stamps/cellranger_vdj/{sample}.stamp")
	log: os.path.abspath("logs/cellranger_vdj/{sample}.log")
	threads: 1
	params: 
		fastqdirs = lambda wildcards: get_count_fastqdirs(wildcards, lib_types={"BCR", "TCR"}, info=info, output_dir=config["output_dir"])
		reference = config["cellranger_vdj_reference"],
		custom_flags = config.get("cellranger_vdj_args", ""),
		output_path = os.path.join(config["output_dir"], "cellranger_vdj")
	envmodules: "cellranger/9.0.0"
	message: "Making BCR/TCR contig file for {wildcards.sample}"
	shell:
		"""
		( \
		mkdir -p stamps/cellranger_vdj && \
		mkdir -p {params.output_path} && \
		cd {params.output_path} && \
		cellranger vdj \
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
cellranger_vdj = [f"stamps/cellranger_vdj/{sample}.stamp" for sample in samples]