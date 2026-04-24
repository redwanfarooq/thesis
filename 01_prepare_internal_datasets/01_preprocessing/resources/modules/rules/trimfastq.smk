##########################################################################################
# Snakemake rule for seqtk trimfq
# Author: Redwan Farooq
# Requires functions from resources/scripts/rule.py
##########################################################################################


# Define rule
rule trimfastq:
	output: os.path.abspath("stamps/trimfastq/{lib}.stamp")
	log: os.path.abspath("logs/trimfastq/{lib}.log")
	threads: 1
	params:
		run_path = lambda wildcards: get_run_path(wildcards, info=info, run_dir=config["run_dir"]),
		lib_type = lambda wildcards: get_lib_type(wildcards, info=info),
		read_trim_flags_R1 = lambda wildcards: get_read_trim_flags(wildcards, read_trim=config.get("read_trim", None), read="R1", info=info),
		read_trim_flags_R2 = lambda wildcards: get_read_trim_flags(wildcards, read_trim=config.get("read_trim", None), read="R2", info=info),
		read_trim_flags_R3 = lambda wildcards: get_read_trim_flags(wildcards, read_trim=config.get("read_trim", None), read="R3", info=info),
		output_path = os.path.join(config["output_dir"], "fastqs") # DO NOT CHANGE - downstream rules will search for FASTQs in this directory
	# conda: "seqtk"
	envmodules: "seqtk"
	message: "Trimming FASTQ files for {wildcards.lib}"
	shell:
		"""
		( \
		mkdir -p stamps/trimfastq && \
		mkdir -p {params.output_path}/{wildcards.lib} && \
		find {params.run_path}/{params.lib_type} -type f -name '*_R1_001.fastq.gz' -exec sh -c 'for file do echo $(date -u +"%F %T") Processing file: $file; seqtk trimfq {params.read_trim_flags_R1} $file | pigz -p {threads} > {params.output_path}/{wildcards.lib}/$(basename $file); done' sh {{}} \; && \
		find {params.run_path}/{params.lib_type} -type f -name '*_R2_001.fastq.gz' -exec sh -c 'for file do echo $(date -u +"%F %T") Processing file: $file; seqtk trimfq {params.read_trim_flags_R2} $file | pigz -p {threads} > {params.output_path}/{wildcards.lib}/$(basename $file); done' sh {{}} \; && \
		find {params.run_path}/{params.lib_type} -type f -name '*_R3_001.fastq.gz' -exec sh -c 'for file do echo $(date -u +"%F %T") Processing file: $file; seqtk trimfq {params.read_trim_flags_R3} $file | pigz -p {threads} > {params.output_path}/{wildcards.lib}/$(basename $file); done' sh {{}} \; && \
		touch {output} \
		) > {log} 2>&1
		"""

# Set rule targets
trimfastq = [f"stamps/trimfastq/{lib}.stamp" for lib in libs if libs[lib]["format"].upper() == "FASTQ" and "read_trim" in config.keys()]