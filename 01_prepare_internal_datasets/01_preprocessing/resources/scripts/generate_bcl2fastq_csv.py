#!/bin/env python


"""
Generates CSV sample sheets for libraries from 10X multiome experiments for use with bcl2fastq.
Requires:
- Metadata table file with the following fields:
    run: run folder name
    format: file type ("BCL" or "FASTQ" - case-insensitive)
    lib_type: library type
    sample_id: sample ID
    sample_index: EITHER index set name OR literal i7 index sequence
    lane: EITHER lane number OR * (for all lanes)
    sample_index2: literal i5 index sequence (if dual indexing used and sample_index is literal i7 index sequence)
Optional:
- Index kit CSV files with the following fields:
    index_name: index set name (e.g. SI-TT-A1); must be first field
    *: 1 or more fields specifying index sequences in the following order:
        Dual index kits: i7 sequence, i5 sequence (forward strand), i5 sequence (reverse complement)
        Single index kits: i7 sequence (additional fields if more than 1 sequence per index set)
"""


# ==============================
# MODULES
# ==============================
import os
import csv
import docopt
from loguru import logger
import pandas as pd
from classes import IndexKit
from id import lib_id


# ==============================
# COMMAND LINE OPTIONS
# ==============================
# Define options
DOC = """
Generate CSV sample sheets for use with bcl2fastq

Usage:
  generate_bcl2fastq_csv.py --md=<md> --outdir=<outdir> [--dual=<dual>] [--single=<single>] [options]

Arguments:
  -m --md=<md>              Metadata table file (required)
  -o --outdir=<outdir>      Output directory (required)
  -d --dual=<dual>          Comma-separated list of dual index kit CSV files
  -s --single=<single>      Comma-separated list of single index kit CSV files

Options:
  -r --reversecomplement    Use reverse complement of i5 index
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
        {"run", "format", "lib_type", "sample_id", "sample_index", "lane"}
    ), "Invalid metadata table file."
    dual = (
        (
            IndexKit(
                file=os.path.basename(file),
                lookup=_read_index_csv(
                    filename=file,
                    # Choose appropriate i5 index sequence for sequencer workflow
                    index_cols=(
                        slice(1, 4, 2) if opt["--reversecomplement"] else slice(1, 3, 1)
                    ),
                ),
                index_type="dual",
            )
            for file in opt["--dual"].split(",")
        )
        if opt["--dual"] is not None
        else None
    )
    single = (
        (
            IndexKit(
                file=os.path.basename(file),
                lookup=_read_index_csv(
                    filename=file,
                    index_cols=slice(1, None, 1),
                ),
                index_type="single",
            )
            for file in opt["--single"].split(",")
        )
        if opt["--single"] is not None
        else None
    )
    kits = (
        {y for x in (dual, single) if x is not None for y in x}
        if any(x is not None for x in (dual, single))
        else None
    )

    # Select runs in BCL format only
    md = md[md.format.str.upper() == "BCL"]

    if not md.empty:
        # Reverse complement literal i5 index sequence if required
        if "sample_index2" in md.columns:
            md["sample_index2"] = md.sample_index2.fillna("")
            if opt["--reversecomplement"]:
                md["sample_index2"] = md.sample_index2.apply(_reverse_complement)

        # Add lane and unique library ID; create a row for each lane if not * and more than one specified
        md = (
            md.assign(
                lane=lambda x: [
                    "" if lane == "*" else str(lane).split() for lane in x.lane
                ],
                lib_id=lambda x: lib_id(x.lib_type.tolist(), x.run.tolist()),
            )
            .explode("lane")
            .reset_index(drop=True)
        )

        # Generate sample sheets
        logger.info("Generating sample sheets for bcl2fastq")
        for x in md.lib_id.unique():
            generate_sample_sheet(
                df=md[md.lib_id == x],
                index_kits=kits,
                filename=os.path.join(opt["--outdir"], f"{x}.csv"),
            )
            logger.success(
                "Output file: {}",
                os.path.abspath(os.path.join(opt["--outdir"], f"{x}.csv")),
            )


def _read_index_csv(filename: str, index_cols: slice) -> dict:
    with open(file=filename, mode="r", encoding="UTF-8") as file:
        data = csv.reader(file)
        return {rows[0]: tuple(rows[index_cols]) for rows in data}


def _get_index_kit(index_names: set[str], index_kits: set[IndexKit]) -> IndexKit | None:
    for index_kit in index_kits:
        if index_kit.match(names=index_names, strict=True):
            return index_kit
    return None


def _reverse_complement(seq: str) -> str:
    if seq != "":
        assert set(seq).issubset("ATCG"), "Invalid bases detected in sequence."
        seq = seq.translate(str.maketrans("ATCG", "TAGC"))[::-1]
    return seq


def generate_sample_sheet(
    df: pd.DataFrame,
    index_kits: set[IndexKit] | None,
    filename: str | None = None,
) -> pd.DataFrame:
    """
    Generate CSV sample sheets for each sequencing run.

    Arguments:
        ``df``: DataFrame containing run metadata for a single library.\n
        ``index_kits``: Set of IndexKit objects to search for matching index kit or ``None``.\n
        ``filename``: Output file path or ``None``.

    Returns:
        Writes formatted CSV string to ``filename`` (if provided) and returns DataFrame containing CSV data.
    """
    index_kit = (
        _get_index_kit(index_names=set(df.sample_index.unique()), index_kits=index_kits)
        if index_kits is not None
        else None
    )

    if index_kit is None:
        out = pd.DataFrame(
            {
                "Lane": df.lane.tolist(),
                "Sample_ID": df.sample_id.tolist(),
                "index": df.sample_index.tolist(),
            }
        )
        if "sample_index2" in df.columns:
            out["index2"] = df.sample_index2.tolist()
    else:
        if index_kit.index_type == "dual":
            out = pd.DataFrame(
                {
                    "Lane": df.lane.tolist(),
                    "Sample_ID": df.sample_id.tolist(),
                    "index": [
                        index_kit.extract(index_name, 0)
                        for index_name in df.sample_index
                    ],
                    "index2": [
                        index_kit.extract(index_name, 1)
                        for index_name in df.sample_index
                    ],
                }
            )
        elif index_kit.index_type == "single":
            out = pd.DataFrame(
                {
                    "Lane": df.lane.repeat(
                        len(index_kit.extract("index_name"))
                    ).tolist(),
                    "Sample_ID": df.sample_id.repeat(
                        len(index_kit.extract("index_name"))
                    ).tolist(),
                    "index": [
                        x
                        for index_name in df.sample_index
                        for x in index_kit.extract(index_name)
                    ],
                }
            )

    if filename is not None:
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        csv_string = "[Data]\n" + out.to_csv(header=True, index=False)
        with open(file=filename, mode="w", encoding="UTF-8") as file:
            file.write(csv_string)

    return out


# ==============================
# SCRIPT
# ==============================
if __name__ == "__main__":
    _main(opt=docopt.docopt(DOC))
