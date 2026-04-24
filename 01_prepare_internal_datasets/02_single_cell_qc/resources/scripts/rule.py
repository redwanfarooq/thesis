"""
Functions for use in Snakemake rule definitions.
"""

import os
import re
import h5py
import yaml


def parse_info(info: dict) -> dict:
    """
    Parse dictionary of sample info.

    Arguments:
        ``info``: dictionary of sample info.

    Returns:
        Dictionary of parsed sample info.
    """
    samples = list(info.keys())

    return {"samples": samples}


def get_merge_flags(wildcards, **kwargs) -> str:
    """
    Get flags for multimodal count matrix merging script.

    Arguments:
        ``wildcards``: Snakemake ``wildcards`` object.
        ``kwargs``: keyword arguments for flags.

    Returns:
        String containing flag to be inserted into shell command.
    """
    flags = [
        f"--{str(key).replace('_', '-')} {str(value).format(sample=wildcards.sample)}"
        for key, value in kwargs.items()
        if value is not None
    ]
    return " ".join(flags)


def get_expected_cells_flag(wildcards, info: dict) -> str:
    """
    Get flag for expected number of cells in CellBender.

    Arguments:
        ``wildcards``: Snakemake ``wildcards`` object.
        ``info``: dictionary of sample info.

    Returns:
        String containing flag to be inserted into shell command.
    """
    n_cells = sum(_.get("cells_loaded", 0) for _ in info[wildcards.sample].values()) * 0.7 # ~70% capture rate of cells loaded (does not need to be precise)
    return f"--expected-cells {round(n_cells)}" if n_cells > 0 else ""


def get_ignore_features_flag(hdf5: str, regex: str | None = None) -> str:
    """
    Get flag for ignoring features in CellBender.

    Arguments:
        ``hdf5``: path to merged multimodal count matrix.
        ``regex``: regular expression to match features to ignore.

    Returns:
        String containing flag to be inserted into shell command.
    """
    if os.path.exists(hdf5) and regex is not None:
        with h5py.File(hdf5, mode="r") as file:
            features = [
                _.decode("UTF-8") for _ in list(file["matrix"]["features"]["id"])
            ]
        features_to_ignore = [
            str(i) for i, _ in enumerate(features) if re.search(regex, _)
        ]
        return f"--ignore-features {' '.join(features_to_ignore)}"
    return ""


def get_features_matrix(
    wildcards, data_dir: str, cellbender: bool = False, filtered: bool | None = False
) -> str:
    """
    Get path to merged multimodal count matrix.

    Arguments:
        ``wildcards``: Snakemake ``wildcards`` object.
        ``path``: path to pipeline data output directory.
        ``cellbender``: boolean indicating whether CellBender is used to preprocess count matrices.
        ``filtered``: boolean indicating whether filtered or raw count matrix is used.

    Returns:
        Path to merged multimodal count matrix.
    """
    return os.path.join(
        data_dir,
        f"{'cellbender' if cellbender else 'merge'}/{wildcards.sample}/{'filtered' if filtered else 'raw'}_feature_bc_matrix.h5",
    )


def get_hto_metadata(wildcards, info: dict) -> str:
    """
    Get cell hashing metadata string.

    Arguments:
        ``wildcards``: Snakemake ``wildcards`` object.
        ``info``: dictionary of sample info.

    Returns:
        YAML-formatted cell hashing metadata string.
    """
    return yaml.dump(info[wildcards.sample])
