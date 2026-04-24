##########################################################################################
# Snakemake rule for count_reads
# Author: Redwan Farooq
# Requires outputs from resources/rules/multiqc.smk
##########################################################################################

# Define rule
rule count_reads:
	input: os.path.abspath("stamps/multiqc/multiqc.stamp")
	output: os.path.abspath("stamps/count_reads/count_reads.stamp")
	log: os.path.abspath("logs/count_reads/count_reads.log")
	threads: 1
	params:
		multiqc_dir = os.path.join(config["output_dir"], "qc/multiqc"),
		output_path = os.path.join(config["output_dir"], "qc/count_reads")
	message: "Calculating total read counts per library type/sample"
	run:
		import pandas as pd
		from loguru import logger
		logger.remove(0)
		logfile = logger.add(sink=log[0], level="INFO", mode="w")

		with logger.catch(reraise=True):
			shell("mkdir -p stamps/count_reads")
			shell("mkdir -p {params.output_path}")

			if os.path.isfile(
				os.path.join(params["multiqc_dir"], "multiqc_data/multiqc_fastqc.txt")
			):
				logger.info("Found source data: {}", os.path.join(params["multiqc_dir"], "multiqc_data/multiqc_fastqc.txt"))
				shell("cp {params.multiqc_dir}/multiqc_data/multiqc_fastqc.txt {params.output_path}")
			elif os.path.isfile(
				os.path.join(params["multiqc_dir"], "multiqc_data.zip")
			):
				logger.info("Found compressed archive: {}", os.path.join(params["multiqc_dir"], "multiqc_data.zip"))
				logger.info("Extracting source data from compressed archive")
				shell("unzip -p {params.multiqc_dir}/multiqc_data.zip multiqc_fastqc.txt > {params.output_path}/multiqc_fastqc.txt")
			else:
				raise FileNotFoundError

			df = pd.read_table(
				os.path.join(params["output_path"], "multiqc_fastqc.txt"),
				delimiter="\t",
			).rename(
				columns={
					"Sample": "directory",
					"Filename": "filename",
					"Total Sequences": "read_count",
				}
			)
			df = (
				df[["directory", "filename", "read_count"]][df.filename.str.contains("_R1_")]
				.assign(
					lib_type=lambda x: [y.split("-")[0] for y in x.directory],
					sample=lambda x: [y.split("_")[0] for y in x.filename],
					read_count=lambda x: [int(y) for y in x.read_count],
				)
				.groupby(["lib_type", "sample"])
				.agg({"read_count": "sum"})
				.reset_index()
				.pivot(index="sample", columns="lib_type", values="read_count")
			)
			df.to_csv(
				os.path.join(params["output_path"], "read_counts.tsv"),
				sep="\t",
				header=True,
				index=True,
			)
			logger.success("Output file: {}", os.path.join(params["output_path"], "read_counts.tsv"))

			shell("rm {params.output_path}/multiqc_fastqc.txt")
			shell("touch {output}")
		
		logger.remove(logfile)


# Set rule targets
count_reads = ["stamps/count_reads/count_reads.stamp"]