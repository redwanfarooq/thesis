##########################################################################################
# Snakemake rule for ambient removal
# Author: Redwan Farooq
# Requires functions from resources/scripts/rule.py
# Requires outputs from resources/rules/droplet_qc.smk
# Requires outputs from resources/rules/libraries_qc.smk
##########################################################################################

scripts_dir = config.get("scripts_dir", "resources/scripts")

# Define rule
rule ambient_removal:
    input: os.path.join(config["output_dir"]["reports"], "libraries_qc/{sample}.html")
    output: os.path.join(config["output_dir"]["reports"], "ambient_removal/{sample}.html")
    log: os.path.abspath("logs/ambient_removal/{sample}.log")
    threads: 1
    params:
        script_path = scripts_dir if os.path.isabs(scripts_dir) else os.path.join(workflow.basedir, scripts_dir),
        droplet_qc = lambda wildcards: os.path.join(config["output_dir"]["data"], "droplet_qc/{sample}"),
        libraries_qc = lambda wildcards: os.path.join(config["output_dir"]["data"], "libraries_qc/{sample}"),
        outdir = lambda wildcards: os.path.join(config["output_dir"]["data"], "finalised/{sample}"),
        features_matrix = lambda wildcards: get_features_matrix(wildcards, data_dir=config["output_dir"]["data"], cellbender="cellbender" in module_rules, filtered=config.get("pre_filtered", None)),
        adt_prefix = config.get("adt_prefix", None),
        gex_ambient_removal = config.get("gex_ambient_removal", None),
        adt_ambient_removal = config.get("adt_ambient_removal", None)
    conda: "single_cell_qc"
    # envmodules: "R-cbrg"
    message: "Running ambient removal for {wildcards.sample}"
    script: "{params.script_path}/ambient_removal.Rmd"

ambient_removal = expand(os.path.join(config["output_dir"]["reports"], "ambient_removal/{sample}.html"), sample=samples)