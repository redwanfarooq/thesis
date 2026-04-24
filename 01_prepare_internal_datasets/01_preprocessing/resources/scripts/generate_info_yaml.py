#!/bin/env python


"""
Generates info YAML for use with preprocessing pipeline.
Requires:
- Metadata table file with the following fields:
    run: run folder name
    format: file type ("BCL" or "FASTQ" - case-insensitive)
    lib_type: library type
    sample_id: sample ID
"""


# ==============================
# MODULES
# ==============================
import os
import yaml
import docopt
from loguru import logger
import pandas as pd
from id import lib_id


# ==============================
# COMMAND LINE OPTIONS
# ==============================
# Define options
DOC = """
Generate info YAML for use with preprocessing pipeline

Usage:
  generate_info_yaml.py --md=<md> --outdir=<outdir> [options]

Arguments:
  -m --md=<md>              Metadata table file (required)
  -o --outdir=<outdir>      Output directory (required)

Options:
  -h --help                 Show this screen
"""


# ==============================
# FUNCTIONS
# ==============================
@logger.catch(reraise=True)
def _main(opt: dict) -> None:
    # Read input CSV and check fields are valid
    md = pd.read_csv(opt["--md"], header=0, sep=None, engine="python")
    assert set(md.columns).issuperset(
        {"run", "format", "lib_type", "sample_id"}
    ), "Invalid metadata table file."

    # Add unique library ID
    md = md.assign(lib_id=lambda x: lib_id(x.lib_type.tolist(), x.run.tolist()))

    # Generate info YAML
    logger.info("Generating info YAML")
    generate_info_yaml(
        df=md[["sample_id", "lib_id", "format", "lib_type", "run"]].drop_duplicates(),
        filename=os.path.join(opt["--outdir"], "info.yaml"),
    )
    logger.success(
        "Output file: {}", os.path.abspath(os.path.join(opt["--outdir"], "info.yaml"))
    )


def generate_info_yaml(df: pd.DataFrame, filename: str | None = None) -> dict:
    """
    Generate info YAML.

    Arguments:
        ``df``: DataFrame containing run metadata.\n
        ``filename``: Output file path or ``None``.

    Returns:
        Writes formatted YAML string to ``filename`` (if provided) and returns dictionary containing YAML data.
    """
    out = {
        sample_id: libs.loc[sample_id].to_dict("index")
        for sample_id, libs in df.set_index(["sample_id", "lib_id"]).groupby(
            "sample_id"
        )
    }

    if filename is not None:
        with open(file=filename, mode="w", encoding="UTF-8") as file:
            yaml.dump(data=out, stream=file)

    return out


# ==============================
# SCRIPT
# ==============================
if __name__ == "__main__":
    _main(opt=docopt.docopt(DOC))
