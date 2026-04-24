"""
Functions for extracting/generating unique IDs for use with preprocessing pipeline scripts.
"""


def paste(*args, sep: str = "") -> list[str]:
    """
    Vectorised string concatenation.

    Arguments:
        ``*args``: 2 or more lists of strings or data type coercible to string.\n
        ``sep``: String used as separator during concatenation. Default: ``""``.

    Returns:
        List of concetenated strings.
    """
    combs = zip(*args)
    return [sep.join(str(j) for j in i) for i in combs]


def lib_id(lib_type: list[str], run: list[str]) -> list[str]:
    """
    Generate unique library IDs from library type and flow cell ID.

    Arguments:
        ``lib_type``: List of strings specifying library types.\n
        ``run``: List of strings specifying Illumina run folder names.

    Returns:
        List of unique library IDs.
    """
    return paste(lib_type, fcid(run), sep="-")


def fcid(run: list[str]) -> list[str]:
    """
    Extract flow cell IDs from Illumina run folder names.

    Arguments:
        ``run``: List of strings specifying Illumina run folder names.

    Returns:
        List of flow cell IDs.
    """
    return map(lambda x: x.split("_")[-1], run)
