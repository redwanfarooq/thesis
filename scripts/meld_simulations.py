#!/bin/env python


"""
Run MELD simulations using randomly permuted sample condition labels to obtain an
empirical null distribution.

Requires precomputed normalised and scaled data matrix, PCA matrix or graph adjacency
matrix and sample metadata table as input.
"""


# ==============================
# MODULES
# ==============================
import gc
import warnings
import gzip
import pathlib
import docopt
import numpy as np
import pandas as pd
import scipy.io
import graphtools as gt
import graphtools.graphs as gt_graphs
import meld
from joblib import Parallel, delayed
from tqdm import tqdm
from loguru import logger
from sklearn.preprocessing import normalize


# ==============================
# COMMAND LINE OPTIONS
# ==============================
# Define options
DOC = """
Run MELD simulations using randomly permuted sample condition labels to obtain an empirical null distribution.

Usage:
    meld_simulations.py (--scaled=<scaled> | --pca=<pca> | --graph=<graph>) --metadata=<metadata> --replicate=<replicate> --condition=<condition> --reference=<reference> --output=<output> [options]

Arguments:
    --scaled=<scaled>               Path to the normalised and scaled expression matrix TSV file (can be compressed)
    --pca=<pca>                     Path to the PCA or embedding matrix TSV file (can be compressed)
    --graph=<graph>                 Path to the graph adjacency matrix in MatrixMarket format
    --metadata=<metadata>           Path to the sample metadata TSV file
    --replicate=<replicate>         Column name for replicate labels in sample metadata
    --condition=<condition>         Column name for condition labels in sample metadata
    --reference=<reference>         Reference condition for computing likelihood ratios
    --output=<output>               Path to the output directory

Options:
    -t --threads=<threads>          Number of threads to use [default: 1]
    -l --log=<log>                  Path to the log file [default: meld_simulations.log]
    --n-pca=<n-pca>                 Number of PCA components to use (ignored if --graph is provided) [default: 50]
    --k=<k>                         Number of neighbors for k-nearest neighbor graph (ignored if --graph is provided) [default: 5]
    --beta=<beta>                   MELD beta parameter [default: 60]
    --n-iter=<n-iter>               Number of simulation iterations [default: 1000]
    -q --quiet                      Do not print progress to stdout
    -h --help                       Show this screen
"""


# ==============================
# FUNCTIONS
# ==============================
@logger.catch(reraise=True)
def _main(opt: dict) -> None:
    # Set up logging
    # Disable logging to stdout if quiet mode
    if opt["--quiet"]:
        logger.remove()
    # Enable logging to file
    logger.add(opt["--log"])
    # Enable logging of warnings
    warnings.showwarning = lambda message, category, filename, lineno, file=None, line=None: logger.opt(
        colors=True
    ).warning(
        "\n\n".join(
            [
                f"{message}",
                f"<n><white>> File '</white><green>{filename}</green><white>', line </white><yellow>{lineno}</yellow></n>",
                f"{category.__name__}<white>: {message}</white>",
            ]
        )
    )

    # Log options
    logger.info("Parsing command line arguments")
    for key in [
        "--scaled",
        "--pca",
        "--graph",
        "--metadata",
        "--replicate",
        "--condition",
        "--reference",
        "--output",
        "--threads",
        "--log",
        "--n-pca",
        "--k",
        "--beta",
        "--n-iter",
        "--quiet",
    ]:
        logger.info(f"{key.replace('--', '')}: {opt[key]}")

    # Load input data
    if opt["--pca"]:
        logger.info(f"Loading PCA matrix: {opt['--pca']}")
        data = pd.read_csv(
            opt["--pca"], sep="\t", header=None, index_col=None
        ).to_numpy(dtype=np.float64)
        data_type = "pca"
        if int(opt["--n-pca"]) <= data.shape[1]:
            data = data[:, : int(opt["--n-pca"])]
            n_pca = None
        else:
            raise ValueError(
                f"Number of PCA components exceeds the number of columns in the input matrix: {data.shape[1]}"
            )
    elif opt["--scaled"]:
        logger.info(f"Loading normalised and scaled data matrix: {opt['--scaled']}")
        data = pd.read_csv(
            opt["--scaled"], sep="\t", header=None, index_col=None
        ).to_numpy(dtype=np.float64)
        data_type = "scaled"
        n_pca = int(opt["--n-pca"])
    elif opt["--graph"]:
        logger.info(f"Loading precomputed graph(s): {opt['--graph']}")
        graph = gt.Graph(
            scipy.io.mmread(opt["--graph"]),
            precomputed="adjacency",
            use_pygsp=True,
        )
        data_type = "graph"
        n_pca = None

    # Load metadata
    logger.info(f"Loading metadata: {opt['--metadata']}")
    metadata = pd.read_csv(
        opt["--metadata"],
        sep="\t",
        index_col=0,
        dtype="category",
    )

    if data_type in {"pca", "scaled"}:
        # Compute k-nearest neighbors graph
        logger.info(f"Computing k-nearest neighbors graph using k = {opt['--k']}")
        graph = gt.Graph(
            data,
            knn=int(opt["--k"]),
            n_pca=n_pca,
            n_jobs=int(opt["--threads"]),
            use_pygsp=True,
            verbose=True,
        )

    # Run simulations in parallel
    logger.info(f"Running {opt['--n-iter']} simulations")
    results = Parallel(n_jobs=int(opt["--threads"]))(
        delayed(_run_simulation)(
            graph=graph,
            beta=int(opt["--beta"]),
            metadata=metadata,
            replicate_col=str(opt["--replicate"]),
            condition_col=str(opt["--condition"]),
            reference=str(opt["--reference"]),
            seed=seed,
        )
        for seed in tqdm(
            range(int(opt["--n-iter"])),
            desc="Progress",
            disable=(int(opt["--threads"]) == 1 or opt["--quiet"]),
        )
    )

    # Aggregate results
    logger.info("Aggregating results")
    results = {
        condition: np.column_stack([result[condition] for result in results])
        for condition in results[0].keys()
    }

    # Save the results
    pathlib.Path(opt["--output"]).mkdir(parents=True, exist_ok=True)
    for condition, data in results.items():
        with gzip.open(
            pathlib.Path(opt["--output"]) / f"{condition}.tsv.gz", "wt"
        ) as f:
            np.savetxt(f, data, delimiter="\t")
    logger.success(f"Output path: {opt['--output']}")


def _permute_condition_labels(
    metadata: pd.DataFrame, replicate_col: str, condition_col: str, seed: int
) -> pd.Series:
    assert isinstance(metadata, pd.DataFrame), "'metadata' must be a DataFrame"
    assert (
        replicate_col in metadata.columns
    ), f"'{replicate_col}' not found in 'metadata'"
    assert (
        condition_col in metadata.columns
    ), f"'{condition_col}' not found in 'metadata'"
    metadata[condition_col] = metadata.groupby(replicate_col)[condition_col].transform(
        lambda x: x.sample(frac=1, random_state=seed).values
    )
    permuted_labels = pd.Series(
        "_".join(x) for x in zip(metadata[replicate_col], metadata[condition_col])
    )
    return permuted_labels


def _compute_log_likelihood_ratio(
    densities: pd.DataFrame, replicate: str, reference: str
) -> pd.DataFrame:
    assert isinstance(densities, pd.DataFrame), "'densities' must be a DataFrame"
    df = densities.filter(regex=f"^{replicate}_")
    normalize(df, norm="l1", axis=1, copy=False)
    reference = df.filter(regex=f"_{reference}$")
    df = np.log(df.div(reference.squeeze(), axis=0))
    df = df.drop(columns=reference.columns)
    return df


def _run_simulation(
    graph: gt_graphs.kNNPyGSPGraph,
    beta: int,
    metadata: pd.DataFrame,
    replicate_col: str,
    condition_col: str,
    reference: str,
    seed: int,
) -> np.ndarray:
    labels = _permute_condition_labels(metadata, replicate_col, condition_col, seed)
    meld_op = meld.MELD(beta=beta, verbose=False)
    densities = meld_op.fit_transform(graph, sample_labels=labels)
    log_lr = [
        _compute_log_likelihood_ratio(densities, replicate, reference)
        for replicate in metadata[replicate_col].unique()
    ]
    log_lr = pd.concat(log_lr, axis=1)
    median_log_lr = {
        condition: log_lr.filter(regex=f"_{condition}$").median(axis=1).values
        for condition in metadata[condition_col].unique()
        if condition != reference
    }
    gc.collect()  # Garbage collection to prevent memory leaks
    return median_log_lr


# ==============================
# SCRIPT
# ==============================
if __name__ == "__main__":
    # Run the main function
    _main(opt=docopt.docopt(DOC))
