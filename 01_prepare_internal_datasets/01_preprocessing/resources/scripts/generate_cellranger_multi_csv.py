#!/bin/env python


"""
Generates CSV configuration sheets for each sample for 10X feature barcoding +/- immune profiling experiments
for use with cellranger multi.
Requires:
- Metadata table file with the following fields:
    run: run folder name
    lib_type: library type
    sample_id: sample ID
- Sample hashing CSV file with the following fields (only if using sample hashing):
    sample_id: sample ID
    hash_id: hash ID
    ocm_barcode_ids/hashtag_ids/cmo_ids: OCM barcode/hashtag/CMO IDs
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
Generate CSV configuration sheets for use with cellranger multi

Usage:
  generate_cellranger_multi_csv.py --md=<md> --fastqdir=<fastqdir> --outdir=<outdir> [--features=<features> --hashes=<hashes> --transcriptome=<transcriptome> --vdj=<vdj>] [options]

Arguments:
  -m --md=<md>                         Metadata table file (required)
  -f --fastqdir=<fastqdir>             FASTQ directory (required)
  -o --outdir=<outdir>                 Output directory (required)
  --features=<features>                Features reference CSV file
  --hashes=<hashes>                    Sample hashing CSV file
  --transcriptome=<transcriptome>      Path to Cell Ranger transcriptome reference
  --vdj=<vdj>                          Path to Cell Ranger VDJ reference
  --create-bam=<bool>                  Enable or disable BAM file generation [default: true]
  --tenx-cloud-token-path=<path>       Path to 10x Cloud Analysis user token used to enable cell annotation
  --cell-annotation-model=<model>      Cell annotation model to use
  --chemistry=<chemistry>              Assay configuration
  --expect-cells=<int>                 Expected number of recovered cells
  --force-cells=<int>                  Force pipeline to use this number of cells
  --include-introns=<bool>             Include intronic reads in gene expression count matrix
  --no-secondary=<bool>                Disable secondary analysis
  --check-library-compatibility=<bool> Evaluate 10x Barcode overlap between libraries when multiple libraries are specified
  --emptydrops-minimum-umis=<int>      Minimum number of UMIs for a cell to be called by EmptyDrops

Options:
  -h --help                            Show this screen
"""


# ==============================
# FUNCTIONS
# ==============================
@logger.catch(reraise=True)
def _main(opt: dict) -> None:
    # Read input CSVs and check fields are valid
    md = pd.read_csv(opt["--md"], header=0, sep=None, engine="python")
    assert set(md.columns).issuperset(
        {"run", "lib_type", "sample_id"}
    ), "Invalid metadata table file."

    if opt["--hashes"]:
        hashes = pd.read_csv(opt["--hashes"], header=0, sep=None, engine="python")
        assert set(hashes.columns).issuperset(
            {"sample_id", "hash_id"}
        ), "Invalid sample hashing CSV file."

    # Add unique library ID
    md = md.assign(lib_id=lambda x: lib_id(x.lib_type.tolist(), x.run.tolist()))

    # Generate library sheets
    logger.info("Generating configuration sheets for cellranger multi")
    for x in md.sample_id.unique():
        generate_config_sheet(
            libraries=md[(md.sample_id == x) & (md.lib_type.isin({"GEX", "ADT", "HTO", "CRISPR", "BCR", "TCR"}))][
                ["sample_id", "lib_id", "lib_type"]
            ].drop_duplicates(),
            fastqdir=opt["--fastqdir"],
            options={k.removeprefix("--"): opt[k] for k in {"--create-bam", "--tenx-cloud-token-path", "--cell-annotation-model", "--chemistry", "--expect-cells", "--force-cells", "--include-introns", "--no-secondary", "--check-library-compatibility", "--emptydrops-minimum-umis"} if k in opt},
            hashes=hashes[hashes.sample_id == x] if opt["--hashes"] else None,
            features=opt["--features"],
            transcriptome=opt["--transcriptome"],
            vdj=opt["--vdj"],
            filename=os.path.join(opt["--outdir"], f"{x}.csv"),
        )
        logger.success(
            "Output file: {}",
            os.path.abspath(os.path.join(opt["--outdir"], f"{x}.csv")),
        )


def _generate_library_sheet(df: pd.DataFrame, fastqdir: str) -> pd.DataFrame:
    out = pd.DataFrame(
        {   
            "fastq_id": df.sample_id.tolist(),
            "fastqs": [os.path.join(fastqdir, lib_id) for lib_id in df.lib_id],
            "feature_types": [
                (
                    "Gene Expression"
                    if lib_type == "GEX"
                    else (
                        "Antibody Capture"
                        if lib_type in {"ADT", "HTO"}
                        else (
                            "CRISPR Guide Capture"
                            if lib_type == "CRISPR"
                            else (
                                "VDJ-B"
                                if lib_type == "BCR"
                                else (
                                    "VDJ-T"
                                    if lib_type == "TCR" else "Unknown"
                                )
                            )
                        )
                    )
                )
                for lib_type in df.lib_type
            ],
        }
    )

    return out

def _generate_sample_sheet(df: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame(
        {
            "sample_id": df.hash_id.tolist(),
            df.columns[-1]: df[df.columns[-1]].tolist(),
        }
    ).groupby("sample_id", as_index=False).agg(lambda x: "|".join(x))

    return out


def generate_config_sheet(
    libraries: pd.DataFrame,
    fastqdir: str,
    options: dict,
    hashes: pd.DataFrame | None = None,
    features: str | None = None,
    transcriptome: str | None = None,
    vdj: str | None = None,
    filename: str | None = None,
) -> str:
    """
    Generate CSV configuration sheets for each sample for 10X feature barcoding +/- immune profiling experiments.

    Arguments:
        ``libraries``: DataFrame containing run metadata for a single sample.\n
        ``fastqdir``: FASTQ directory.\n
        ``hashes``: DataFrame containing sample hashing metadata for a single sample.\n
        ``features``: Features reference CSV file path or ``None``.\n
        ``transcriptome``: Cell Ranger transcriptome reference path or ``None``.\n
        ``vdj``: Cell Ranger VDJ reference path or ``None``.\n
        ``filename``: Output file path or ``None``.

    Returns:
        Writes formatted CSV string to ``filename`` (if provided) and returns CSV string.
    """
    libraries = _generate_library_sheet(libraries, fastqdir)
    if hashes is not None:
        hashes = _generate_sample_sheet(hashes) if not hashes.empty else None

    out = ["[libraries]"]
    out.append(libraries.to_csv(header=True, index=False))
    if hashes is not None:
        out.append("[samples]")
        out.append(hashes.to_csv(header=True, index=False))
    if any(libraries.feature_types == "Gene Expression"):
        assert transcriptome, "Transcriptome reference must be provided for gene expression libraries."
        out.append("[gene-expression]")
        out.append(f"reference,{transcriptome}")
        for k, v in options.items():
            if v is not None:
                out.append(f"{k},{v}")
        out.append("")
    if any(libraries.feature_types.isin({"VDJ-B", "VDJ-T"})):
        assert vdj, "VDJ reference must be provided for VDJ libraries."
        out.append("[vdj]")
        out.append(f"reference,{vdj}\n")
    if any(libraries.feature_types.isin({"Antibody Capture", "CRISPR Guide Capture"})):
        assert features, "Features reference must be provided for antibody capture or CRISPR guide capture libraries."
        out.append("[feature]")
        out.append(f"reference,{features}\n")
    out = "\n".join(out)

    if filename:
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        with open(file=filename, mode="w", encoding="UTF-8") as file:
            file.write(out)

    return out

# ==============================
# SCRIPT
# ==============================
if __name__ == "__main__":
    _main(opt=docopt.docopt(DOC))
