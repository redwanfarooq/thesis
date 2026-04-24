##########################################################################################
# Snakemake rule for creating symbolic links to FASTQ files in the output directory
# Author: Redwan Farooq
# Requires functions from resources/scripts/rule.py
##########################################################################################


# Define rule
rule linkfastq:
	output: os.path.abspath("stamps/linkfastq/{lib}.stamp")
	log: os.path.abspath("logs/linkfastq/{lib}.log")
	threads: 1
	params:
		run_path = lambda wildcards: get_run_path(wildcards, info=info, run_dir=config["run_dir"]),
		lib_type = lambda wildcards: get_lib_type(wildcards, info=info),
		output_path = os.path.join(config["output_dir"], "fastqs") # DO NOT CHANGE - downstream rules will search for FASTQs in this directory
	message: "Creating symbolic links to FASTQ files for {wildcards.lib}"
	shell:
		"""
		( \
		mkdir -p stamps/linkfastq && \
		mkdir -p {params.output_path}/{wildcards.lib} && \
		find {params.run_path}/{params.lib_type} -type f -name '*_R*_001.fastq.gz' -exec sh -c 'for file do echo $(date -u +"%F %T") Processing file: $file; ln -s $file {params.output_path}/{wildcards.lib}/$(basename $file); done' sh {{}} \; && \
		touch {output} \
		) > {log} 2>&1
		"""

# Set rule targets
linkfastq = [f"stamps/linkfastq/{lib}.stamp" for lib in libs if libs[lib]["format"].upper() == "FASTQ" and "read_trim" not in config.keys()]