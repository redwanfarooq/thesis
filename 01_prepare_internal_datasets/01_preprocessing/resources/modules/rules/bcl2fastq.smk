##########################################################################################
# Snakemake rule for bcl2fastq
# Author: Redwan Farooq
# Requires functions from resources/scripts/rule.py
##########################################################################################

# Define rule
rule bcl2fastq:
	output: os.path.abspath("stamps/bcl2fastq/{lib}.stamp")
	log: os.path.abspath("logs/bcl2fastq/{lib}.log")
	threads: 1
	params:
		run_path = lambda wildcards: get_run_path(wildcards, info=info, run_dir=config["run_dir"]),
		samplesheet_path = os.path.abspath(os.path.join(config.get("metadata_dir", "metadata"), "bcl2fastq")),
		bases_mask_flag = lambda wildcards: get_bases_mask_flag(wildcards, bases_mask=config.get("bases_mask", None), info=info, run_dir=config["run_dir"]),
		custom_flags = config.get("bcl2fastq_args", ""),
		output_path = os.path.join(config["output_dir"], "fastqs") # DO NOT CHANGE - downstream rules will search for FASTQs in this directory
	# conda: "bcl2fastq"
	envmodules: "bcl2fastq/2.20.0.422"
	message: "Making FASTQ files for {wildcards.lib}"
	shell:
		"""
		( \
		mkdir -p stamps/bcl2fastq && \
		mkdir -p {params.output_path}/{wildcards.lib} && \
		bcl2fastq \
			--runfolder-dir={params.run_path} \
			--sample-sheet={params.samplesheet_path}/{wildcards.lib}.csv \
			--processing-threads={threads} \
			{params.bases_mask_flag} {params.custom_flags} \
			--output-dir={params.output_path}/{wildcards.lib} && \
		find {params.output_path}/{wildcards.lib} -type f -name 'Undetermined_S0_*.fastq.gz' -exec rm -rf {{}} \; && \
		touch {output} \
		) > {log} 2>&1
		"""

	
# Set rule targets
bcl2fastq = [f"stamps/bcl2fastq/{lib}.stamp" for lib in libs if libs[lib]["format"].upper() == "BCL"]