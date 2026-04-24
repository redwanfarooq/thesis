#!/bin/env python


"""
Generates CSV library sheets (GEX and ATAC) for each sample from 10X multiome experiments
for use with cellranger-arc count.
Requires:
- Metadata table file with the following fields:
    run: run folder name
    lib_type: library type
    sample_id: sample ID
"""


# ==============================
# MODULES
# ==============================
import os
import docopt
from loguru import logger
import pandas as pd
from id import lib_id


# ==============================
# COMMAND LINE OPTIONS
# ==============================
# Define options
DOC = """
Generate CSV GEX and ATAC library sheets for use with cellranger-arc count

Usage:
  generate_cellranger_arc_csv.py --md=<md> --fastqdir=<fastqdir> --outdir=<outdir> [options]

Arguments:
  -m --md=<md>              Metadata table file (required)
  -f --fastqdir=<fastqdir>  FASTQ directory (required)
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
        {"run", "lib_type", "sample_id"}
    ), "Invalid metadata table file."

    # Add unique library ID
    md = md.assign(lib_id=lambda x: lib_id(x.lib_type.tolist(), x.run.tolist()))

    # Generate library sheets
    logger.info("Generating library sheets for cellranger-arc count")
    for x in md.sample_id.unique():
        generate_library_sheet(
            df=md[(md.sample_id == x) & (md.lib_type.isin({"GEX", "ATAC"}))][
                ["sample_id", "lib_id", "lib_type"]
            ].drop_duplicates(),
            fastqdir=opt["--fastqdir"],
            filename=os.path.join(opt["--outdir"], f"{x}.csv"),
        )
        logger.success(
            "Output file: {}",
            os.path.abspath(os.path.join(opt["--outdir"], f"{x}.csv")),
        )


def generate_library_sheet(
    df: pd.DataFrame, fastqdir: str, filename: str | None = None
) -> pd.DataFrame:
    """
    Generate CSV library sheets (GEX and ATAC) for each sample from a 10X multiome experiment.

    Arguments:
        ``df``: DataFrame containing run metadata for a single sample.\n
        ``fastqdir``: FASTQ directory.\n
        ``filename``: Output file path or ``None``.

    Returns:
        Writes formatted CSV string to ``filename`` (if provided) and returns DataFrame containing CSV data.
    """
    out = pd.DataFrame(
        {
            "fastqs": [os.path.join(fastqdir, lib_id) for lib_id in df.lib_id],
            "sample": df.sample_id.tolist(),
            "library_type": [
                "Gene Expression" if lib_type == "GEX" else "Chromatin Accessibility"
                for lib_type in df.lib_type
            ],
        }
    )

    if filename is not None:
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        with open(file=filename, mode="w", encoding="UTF-8") as file:
            out.to_csv(path_or_buf=file, header=True, index=False)

    return out


# ==============================
# SCRIPT
# ==============================
if __name__ == "__main__":
    _main(opt=docopt.docopt(DOC))
