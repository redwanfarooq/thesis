"""
Wrapper function for COMPOSITE multiplet calling algorithm.
"""

import pandas as pd
import sccomposite


def run_composite(file: str, **kwargs) -> None:
    """
    Run COMPOSITE multiplet calling algorithm.

    Arguments:
        file: str: Output file path.
        **kwargs: dict: Key-value pairs of modality and data matrix.

    Returns:
        Writes the output to a file in tab-separated format and returns None.
    """
    if len(kwargs.keys()) == 1:
        if "RNA" in kwargs:
            x, y = sccomposite.RNA_modality.composite_rna(kwargs["RNA"])
        elif "ATAC" in kwargs:
            x, y = sccomposite.ATAC_modality.composite_atac(kwargs["ATAC"])
        elif "ADT" in kwargs:
            x, y = sccomposite.ADT_modality.composite_adt(kwargs["ADT"])
        else:
            raise ValueError("Invalid modality.")
    else:
        x, y = sccomposite.Multiomics.composite_multiomics(**kwargs)

    return pd.DataFrame({"x": x, "y": y}).to_csv(
        file, sep="\t", index=False, header=False
    )
