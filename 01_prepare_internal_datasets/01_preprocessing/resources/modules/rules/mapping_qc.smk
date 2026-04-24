##########################################################################################
# Snakemake rule for mapping QC
# Author: Redwan Farooq
# Requires outputs from resources/rules/starsolo.smk
# Requires outputs from resources/rules/barcounter.smk
# Requires outputs from resources/rules/chromap.smk
# Requires outputs from resources/rules/macs2.smk
##########################################################################################

scripts_dir = config.get("scripts_dir", "resources/scripts")

# Define rule
rule mapping_qc:
	input: [os.path.abspath(path) for path in expand("stamps/{rule}/{sample}.stamp", rule=[rule for rule in module_rules if rule in {"starsolo", "barcounter", "chromap", "macs2"}], sample=samples)]
	output: os.path.abspath("stamps/mapping_qc/mapping_qc.stamp")
	log: os.path.abspath("logs/mapping_qc/mapping_qc.log")
	threads: 1
	params:
		script_path = scripts_dir if os.path.isabs(scripts_dir) else os.path.join(workflow.basedir, scripts_dir),
		input_dir = config["output_dir"],
		samples = ",".join(samples),
		output_path = os.path.join(config["output_dir"], "qc/mapping_qc")
	conda: "quarto"
	# envmodules: "python-cbrg"
	message: "Generating mapping QC report"
	shell:
		"""
		( \
		mkdir -p stamps/mapping_qc && \
		mkdir -p {params.output_path} && \
		cp {params.script_path}/mapping_qc_report.qmd {params.output_path} && \
		cd {params.output_path} && \
		quarto render mapping_qc_report.qmd \
			-P input_dir:{params.input_dir} \
			-P samples:{params.samples} && \
		rm mapping_qc_report.qmd && \
		touch {output} \
		) > {log} 2>&1
		"""


# Set rule targets
mapping_qc = ["stamps/mapping_qc/mapping_qc.stamp"]