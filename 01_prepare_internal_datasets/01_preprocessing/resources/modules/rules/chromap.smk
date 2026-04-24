##########################################################################################
# Snakemake rule for chromap
# Author: Redwan Farooq
# Requires functions from resources/scripts/rule.py
# Requires outputs from resources/rules/bcl2fastq.smk
##########################################################################################

# Define rule
rule chromap:
	input: lambda wildcards: get_count_inputs(wildcards, lib_types={"ATAC"}, info=info, read_trim=True if "read_trim" in config.keys() else False),
	output: os.path.abspath("stamps/chromap/{sample}.stamp")
	log: os.path.abspath("logs/chromap/{sample}.log")
	threads: 1
	params: 
		R1_fastqs = lambda wildcards: get_count_fastqs(wildcards, lib_types={"ATAC"}, read="R1", info=info, output_dir=config["output_dir"]),
		R2_fastqs = lambda wildcards: get_count_fastqs(wildcards, lib_types={"ATAC"}, read="R2", info=info, output_dir=config["output_dir"]),
		R3_fastqs = lambda wildcards: get_count_fastqs(wildcards, lib_types={"ATAC"}, read="R3", info=info, output_dir=config["output_dir"]),
		index = config["chromap_index"],
		reference = config["chromap_reference"],
		whitelist = config["atac_barcode_whitelist"],
		custom_flags = config.get("chromap_args", ""),
		output_path = os.path.join(config["output_dir"], "chromap_macs2") # DO NOT CHANGE - downstream rules will search for fragment files and mapping statistics in this directory
	conda: "chromap"
	# envmodules:
	# 	"chromap/0.2.5",
	# 	"htslib/1.18"
	message: "Making ATAC fragments file for {wildcards.sample}"
	shell:
		"""
		( \
		mkdir -p stamps/chromap && \
		mkdir -p {params.output_path}/{wildcards.sample} && \
		chromap \
			-1 {params.R1_fastqs} \
			-2 {params.R3_fastqs} \
			-b {params.R2_fastqs} \
			-x {params.index} \
			-r {params.reference} \
			--barcode-whitelist <(zcat {params.whitelist}) \
			-o {params.output_path}/{wildcards.sample}/fragments.tsv \
			--summary {params.output_path}/{wildcards.sample}/chromap_summary.csv \
			{params.custom_flags} \
			-t {threads} && \
		bgzip -f -@ {threads} {params.output_path}/{wildcards.sample}/fragments.tsv && \
		tabix -p bed {params.output_path}/{wildcards.sample}/fragments.tsv.gz && \
		touch {output} \
		) > {log} 2>&1
		cp {log} {params.output_path}/{wildcards.sample}/chromap.out
		"""


# Set rule targets
chromap = [f"stamps/chromap/{sample}.stamp" for sample in samples]