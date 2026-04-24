##########################################################################################
# Snakemake rule for summary
# Author: Redwan Farooq
# Requires outputs from resources/rules/droplet_qc.smk
# Requires outputs from resources/rules/libraries_qc.smk
##########################################################################################

scripts_dir = config.get("scripts_dir", "resources/scripts")

# Define rule
rule summary:
    input: expand(os.path.join(config["output_dir"]["reports"], "libraries_qc/{sample}.html"), sample=samples)
    output: os.path.join(config["output_dir"]["reports"], "summary/summary_report.html")
    log: os.path.abspath("logs/summary/summary.log")
    threads: 1
    params:
        script_path = scripts_dir if os.path.isabs(scripts_dir) else os.path.join(workflow.basedir, scripts_dir),
        input_dir = config["output_dir"]["data"],
        samples = ",".join(samples),
        read_counts = config.get("read_counts", ""),
        target_depth = ",".join([f"{k}:{v}" for k, v in config.get("target_depth", {"GEX": 20000, "ATAC": 25000, "ADT": 5000, "HTO": 2000, "CRISPR": 5000, "BCR": 5000, "TCR": 5000}).items()]),
        output_path = os.path.join(config["output_dir"]["reports"], "summary")
    conda: "quarto"
    # envmodules: "python-cbrg"
    message: "Generating summary report"
    shell: 
        """
        ( \
		mkdir -p {params.output_path} && \
		cp {params.script_path}/summary_report.qmd {params.output_path} && \
		cd {params.output_path} && \
		quarto render summary_report.qmd \
			-P input_dir:{params.input_dir} \
			-P samples:{params.samples} \
            -P read_counts:{params.read_counts} \
            -P target_depth:{params.target_depth} && \
		rm summary_report.qmd \
		) > {log} 2>&1
        """

summary = [os.path.join(config["output_dir"]["reports"], "summary/summary_report.html")]