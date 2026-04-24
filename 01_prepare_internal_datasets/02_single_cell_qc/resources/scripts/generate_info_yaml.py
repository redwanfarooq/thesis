#!/bin/env python


"""
Generates info YAML for use with single cell data QC pipeline.
Requires:
- Metadata table file with the following fields:
    sample_id: unique sample ID
    hto: hashtag antibody ID (or 'None' if not used)
    cells_loaded: number of cells loaded
    hash_id: hash ID (or 'None' if not used)
"""


# ==============================
# MODULES
# ==============================
import os
import yaml
import docopt
from loguru import logger
import pandas as pd


# ==============================
# COMMAND LINE OPTIONS
# ==============================
# Define options
DOC = """
Generate info YAML for use with single cell data QC pipeline

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
        {"sample_id", "hto", "cells_loaded", "hash_id"}
    ), "Invalid metadata table file."

    # Generate info YAML
    logger.info("Generating info YAML")
    generate_info_yaml(
        df=md[["sample_id", "hto", "cells_loaded", "hash_id"]],
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
        sample_id: hto.loc[sample_id].to_dict("index")
        for sample_id, hto in df.set_index(["sample_id", "hto"]).groupby("sample_id")
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
