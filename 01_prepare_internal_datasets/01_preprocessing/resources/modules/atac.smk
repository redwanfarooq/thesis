##########################################################################################
# Snakemake module
# Author: Redwan Farooq
# Module name: atac
##########################################################################################


# --------------------------------------------------
# SETUP
# Load modules
import os
import yaml
from resources.scripts.rule import *

# Load and parse sample/library info from YAML file
with open(file=os.path.join(config.get("metadata_dir", "metadata"), "info.yaml"), mode="r", encoding="UTF-8") as file:
    info = yaml.load(stream=file, Loader=yaml.SafeLoader)
for key, value in parse_info(info).items():
    globals()[key] = value

# Set module rules list
module_rules = ['bcl2fastq', 'trimfastq', 'linkfastq', 'fastqc', 'multiqc', 'count_reads', 'chromap', 'macs2', 'mapping_qc']

# Import rules
include: 'rules/bcl2fastq.smk'
include: 'rules/trimfastq.smk'
include: 'rules/linkfastq.smk'
include: 'rules/fastqc.smk'
include: 'rules/multiqc.smk'
include: 'rules/count_reads.smk'
include: 'rules/chromap.smk'
include: 'rules/macs2.smk'
include: 'rules/mapping_qc.smk'

# Set targets list
targets = [x for rule in [bcl2fastq, trimfastq, linkfastq, fastqc, multiqc, count_reads, chromap, macs2, mapping_qc] for x in rule]
# --------------------------------------------------


# --------------------------------------------------
# RULES
rule all:
	input: [os.path.abspath(x) for x in targets]
# --------------------------------------------------