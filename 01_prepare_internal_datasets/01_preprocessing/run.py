#!/bin/env python


"""
Runs preprocessing pipeline.
"""


# ==============================
# MODULES
# ==============================
import os
import glob
import sys
import shutil
import hashlib
from datetime import datetime
import yaml
import docopt
import pandas as pd
from loguru import logger


# ==============================
# COMMAND LINE OPTIONS
# ==============================
# Define options
DOC = """
Run preprocessing pipeline

Usage:
  run.py [options]

Options:
  -h --help                 Show this screen
  -c --clean                Remove existing pipeline results in output directory
  -f --force                Overwrite existing pipeline results in output directory
  -u --update               Update module scripts using rule specifications in 'config/modules.yaml' (will not run pipeline)
"""


# ==============================
# FUNCTIONS
# ==============================
@logger.catch(reraise=True)
def _main(opt: dict) -> None:
    # Get and execute shell command
    cmd = _get_cmd(update=opt["--update"])
    if not opt["--update"]:
        # Check if output directory contains results of a previous pipeline run
        if os.path.exists(os.path.join(OUTPUT_DIR, ".pipeline", "md5sum")):
            if _file_to_str(
                os.path.join(OUTPUT_DIR, ".pipeline", "md5sum")
            ) != _get_hash(config, "VERSION", *METADATA):
                if opt["--clean"]:
                    logger.critical(
                        "Removing existing pipeline results in {}", OUTPUT_DIR
                    )
                    match str(
                        input("Are you sure you want to continue? (y/N) ")
                    ).lower():
                        case "y" | "yes":
                            shutil.rmtree(OUTPUT_DIR)
                        case _:
                            logger.error("Pipeline aborted")
                            sys.exit(1)
                elif opt["--force"]:
                    logger.critical(
                        "Overwriting existing pipeline results in {}", OUTPUT_DIR
                    )
                    match str(
                        input("Are you sure you want to continue? (y/N) ")
                    ).lower():
                        case "y" | "yes":
                            pass
                        case _:
                            logger.error("Pipeline aborted")
                            sys.exit(1)
                else:
                    logger.warning(
                        "Specified output directory contains results of a previous pipeline run from {} with different version, configuration and/or metadata ({}).\nUse --force to overwrite existing results (if required).\nUse --clean to remove existing results.",
                        _file_to_str(
                            os.path.join(OUTPUT_DIR, ".pipeline", "timestamp")
                        ),
                        os.path.join(OUTPUT_DIR, ".pipeline"),
                    )
                    sys.exit(0)
        logger.info("Starting pipeline using module {}", MODULE)
    if not os.system(" && ".join(cmd)) and not opt["--update"]:
        # Save run configuration and metadata to output directory for reproducibility
        os.makedirs(os.path.join(OUTPUT_DIR, ".pipeline"), exist_ok=True)
        with open(
            file=os.path.join(OUTPUT_DIR, ".pipeline", "md5sum"),
            mode="w",
            encoding="UTF-8",
        ) as f:
            f.write(_get_hash(config, "VERSION", *METADATA))
        with open(
            file=os.path.join(OUTPUT_DIR, ".pipeline", "timestamp"),
            mode="w",
            encoding="UTF-8",
        ) as f:
            f.write(str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        shutil.copy(src=CONFIG, dst=os.path.join(OUTPUT_DIR, ".pipeline"))
        shutil.copy(src="VERSION", dst=os.path.join(OUTPUT_DIR, ".pipeline"))
        for f in METADATA:
            shutil.copy(src=f, dst=os.path.join(OUTPUT_DIR, ".pipeline"))
        if os.path.exists("logs"):
            shutil.copytree(
                src="logs",
                dst=os.path.join(OUTPUT_DIR, ".pipeline", "logs"),
                dirs_exist_ok=True,
            )
        # Clean up cluster logs
        if glob.glob("sps-*"):
            os.system("rm -r sps-*")
        logger.success(
            "Pipeline completed successfully; run information saved to {}",
            os.path.join(OUTPUT_DIR, ".pipeline"),
        )


def _cmd(*args):
    cmd = [" ".join(args)]
    return cmd


def _get_index_kit_flags(dual: list[str] | None, single: list[str] | None) -> str:
    dual_flag = f"--dual={','.join(dual)}" if dual is not None else ""
    single_flag = f"--single={','.join(single)}" if single is not None else ""
    return f"{dual_flag} {single_flag}"


def _get_reverse_complement_flag(reverse_complement: bool) -> str:
    return "--reversecomplement" if reverse_complement else ""


def _get_options_flags(options: dict) -> str:
    return " ".join(
        f"--{k}={v if isinstance(v, str) else str(v).lower()}"
        for k, v in options.items()
        if options
    )


def _get_cmd(update: bool = False) -> list[str]:
    if update:
        cmd = _cmd(
            f"{SCRIPTS_DIR}/generate_modules.py",
            "--modules=config/modules.yaml",
            "--template=resources/templates/module.template",
            "--outdir=resources/modules",
        )
    else:
        cmd = _cmd(
            f"{SCRIPTS_DIR}/generate_wrapper.py",
            f"--module={MODULE}",
            "--template=resources/templates/wrapper.template",
        )
        if "bcl2fastq" in RULES and "BCL" in INPUT_FORMATS:
            cmd += _cmd(
                f"{SCRIPTS_DIR}/generate_bcl2fastq_csv.py",
                f"--md={INPUT_TABLE}",
                f"--outdir={os.path.join(METADATA_DIR, 'bcl2fastq')}",
                _get_index_kit_flags(
                    dual=DUAL,
                    single=SINGLE,
                ),
                _get_reverse_complement_flag(reverse_complement=REVERSE_COMPLEMENT),
            )
        if "cellranger" in RULES:
            cmd += _cmd(
                f"{SCRIPTS_DIR}/generate_cellranger_csv.py",
                f"--md={INPUT_TABLE}",
                f"--fastqdir={os.path.join(OUTPUT_DIR, 'fastqs')}",
                f"--outdir={os.path.join(METADATA_DIR, 'cellranger')}",
            )
        if "cellranger_arc" in RULES:
            cmd += _cmd(
                f"{SCRIPTS_DIR}/generate_cellranger_arc_csv.py",
                f"--md={INPUT_TABLE}",
                f"--fastqdir={os.path.join(OUTPUT_DIR, 'fastqs')}",
                f"--outdir={os.path.join(METADATA_DIR, 'cellranger_arc')}",
            )
        if "cellranger_multi" in RULES:
            cmd += _cmd(
                f"{SCRIPTS_DIR}/generate_cellranger_multi_csv.py",
                f"--md={INPUT_TABLE}",
                f"--fastqdir={os.path.join(OUTPUT_DIR, 'fastqs')}",
                f"--outdir={os.path.join(METADATA_DIR, 'cellranger')}",
                f"--features={os.path.abspath(FEATURES)}" if FEATURES else "",
                f"--hashes={os.path.abspath(HASHES)}" if HASHES else "",
                f"--transcriptome={TRANSCRIPTOME}" if TRANSCRIPTOME else "",
                f"--vdj={VDJ}" if VDJ else "",
                _get_options_flags(options=CELLRANGER_MULTI_OPTIONS),
            )
        cmd += _cmd(
            f"{SCRIPTS_DIR}/generate_info_yaml.py",
            f"--md={INPUT_TABLE}",
            f"--outdir={METADATA_DIR}",
        )
        cmd += _cmd("snakemake --profile=profile")
    return cmd


def _get_hash(options: dict, *args):
    x = [_file_to_str(_) for _ in args] if args else []
    x.append(yaml.dump(options, sort_keys=True))
    return hashlib.md5("".join(x).encode()).hexdigest()


def _file_to_str(path: os.PathLike) -> str:
    if os.path.isfile(path):
        with open(file=path, mode="r", encoding="UTF-8") as f:
            return f.read().strip()
    return ""


# ==============================
# SCRIPT
# ==============================
CONFIG = "config/config.yaml"
with open(file=CONFIG, mode="r", encoding="UTF-8") as file:
    config = yaml.load(stream=file, Loader=yaml.SafeLoader)
    SCRIPTS_DIR = config.get("scripts_dir", "resources/scripts")
    METADATA_DIR = config.get("metadata_dir", "metadata")
    TAGS = (
        os.path.join(METADATA_DIR, config["tags"])
        if config.get("tags", None)
        else None
    )
    FEATURES = (
        os.path.join(METADATA_DIR, config["features"])
        if config.get("features", None)
        else None
    )
    HASHES = (
        os.path.join(METADATA_DIR, config["hashes"])
        if config.get("hashes", None)
        else None
    )
    DUAL = config.get("dual_index_kits", None)
    SINGLE = config.get("single_index_kits", None)
    REVERSE_COMPLEMENT = config.get("reverse_complement", False)
    TRANSCRIPTOME = config.get("cellranger_reference", None)
    VDJ = config.get("cellranger_vdj_reference", None)
    CELLRANGER_MULTI_OPTIONS = config.get("cellranger_multi_options", {})
    try:
        INPUT_TABLE = os.path.join(METADATA_DIR, config["runs"])
        df = pd.read_csv(INPUT_TABLE, header=0, sep=None, engine="python")
        assert "format" in df.columns
        INPUT_FORMATS = set(df.format.dropna().str.upper().unique())
        OUTPUT_DIR = config["output_dir"]
        MODULE = config["module"]
    except KeyError as err:
        logger.exception("{} not specified in {}", err, file.name)
        raise KeyError from err
    except FileNotFoundError as err:
        logger.exception("Input table not found: {}", INPUT_TABLE)
        raise FileNotFoundError from err
    except AssertionError as err:
        logger.exception("Input table does not contain 'format' column: {}", INPUT_TABLE)
        raise AssertionError from err
METADATA = [
    _ for _ in [INPUT_TABLE, TAGS, FEATURES, HASHES] if _ is not None and os.path.isfile(_)
]

with open(file="config/modules.yaml", mode="r", encoding="UTF-8") as file:
    try:
        RULES = yaml.load(stream=file, Loader=yaml.SafeLoader)[MODULE]
    except KeyError as err:
        logger.exception("Module {} not specified in {}", err, file.name)
        raise KeyError from err

if __name__ == "__main__":
    _main(opt=docopt.docopt(DOC))
