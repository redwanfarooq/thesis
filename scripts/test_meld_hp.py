#!/bin/env python


"""
Find optimal MELD hyperparameters using grid search and simulations to test performance
of parameter combinations.

Requires precomputed normalised and scaled data matrix, PCA matrix or graph adjacency
matrices as input.
"""


# ==============================
# MODULES
# ==============================
import gc
import warnings
import pathlib
import itertools
import numpy as np
import pandas as pd
import scipy
import graphtools as gt
import graphtools.graphs as gt_graphs
import docopt
import meld
from loguru import logger
from joblib import Parallel, delayed
from tqdm import tqdm


# ==============================
# COMMAND LINE OPTIONS
# ==============================
# Define options
DOC = """
Find optimal MELD hyperparameters using grid search and simulations to test performance
of parameter combinations.

Requires precomputed normalised and scaled data matrix, PCA matrix or graph adjacency
matrices as input.

Usage:
    test_meld_hp.py (--scaled=<matrix> | --pca=<matrix> | --graph=<matrix>) --output=<output> [options]

Arguments:
    --scaled=<matrix>               Path to the normalised and scaled data matrix TSV file (can be compressed)
    --pca=<matrix>                  Path to the PCA matrix TSV file (can be compressed)
    --graph=<matrix>                Semicolon-separated paths to the graph adjacency matrices in MatrixMarket format
    --output=<output>               Path to the output TSV file

Options:
    -t --threads=<threads>          Number of threads to use [default: 1]
    -l --log=<log>                  Path to the log file [default: test_meld_hp.log]
    --seed=<seed>                   Random seed [default: 42]
    --n-pca=<n-pca>                 Number of PCA components to use [default: 50]
    --k-lower=<k-lower>             Lower bound 'k' for k-nearest neighbors graph [default: 1]
    --k-upper=<k-upper>             Upper bound 'k' for k-nearest neighbors graph [default: 26]
    --k-step=<k-step>               Step size 'k' for k-nearest neighbors graph [default: 1]
    --b-lower=<b-lower>             Lower bound 'beta' parameter [default: 1]
    --b-upper=<b-upper>             Upper bound 'beta' parameter [default: 201]
    --b-step=<b-step>               Step size 'beta' parameter [default: 1]
    --n-iter=<n-iter>               Number of simulation iterations for each parameter combination [default: 25]
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
        "--output",
        "--threads",
        "--log",
        "--seed",
        "--n-pca",
        "--k-lower",
        "--k-upper",
        "--k-step",
        "--b-lower",
        "--b-upper",
        "--b-step",
        "--n-iter",
        "--quiet",
    ]:
        logger.info(f"{key.replace('--', '')}: {opt[key]}")

    # Define the search space
    k_range = np.arange(
        int(opt["--k-lower"]), int(opt["--k-upper"]), int(opt["--k-step"])
    )
    beta_range = np.arange(
        int(opt["--b-lower"]), int(opt["--b-upper"]), int(opt["--b-step"])
    )

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
        graphs = tuple(
            gt.Graph(
                scipy.io.mmread(file),
                precomputed="adjacency",
                random_state=int(opt["--seed"]),
                use_pygsp=True,
            )
            for file in opt["--graph"].split(";")
        )
        assert len(k_range) == len(
            graphs
        ), "Number of graphs does not match the number of k values"
        data = dict(zip(k_range, graphs))
        data_type = "graph"
        n_pca = None

    # Run the grid search
    results = grid_search(
        data=data,
        data_type=data_type,
        n_pca=n_pca,
        k_range=k_range,
        beta_range=beta_range,
        n_iter=int(opt["--n-iter"]),
        seed=int(opt["--seed"]),
        n_threads=int(opt["--threads"]),
        progress=not opt["--quiet"],
    )

    # Save the results
    pathlib.Path(opt["--output"]).parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(opt["--output"], sep="\t", index=False)
    logger.success(f"Output file: {opt['--output']}")


def _compute_graph(
    data: np.ndarray, n_pca: int | None, k: int, seed: int
) -> tuple[int, gt_graphs.kNNPyGSPGraph]:
    graph = gt.Graph(
        data=data,
        n_pca=n_pca,
        knn=k,
        random_state=seed,
        use_pygsp=True,
    )
    gc.collect()  # garbage collection to prevent memory leaks
    return k, graph


def _initialise_benchmarker(
    graph: gt_graphs.kNNPyGSPGraph, k: int, n_threads: int, seed: int
) -> tuple[int, meld.Benchmarker]:
    bm = meld.Benchmarker()
    bm.fit_phate(graph, n_jobs=n_threads, random_state=seed, verbose=False)
    bm.graph = graph
    return k, bm


def _run_simulation(
    k: int,
    bm: meld.Benchmarker,
    beta: int,
    iter: int,
) -> tuple[float, int, int, int]:
    bm.set_seed(iter)
    bm.generate_ground_truth_pdf()
    bm.generate_sample_labels()
    bm.calculate_MELD_likelihood(beta=beta)
    mse = bm.calculate_mse(bm.expt_likelihood)
    gc.collect()  # garbage collection to prevent memory leaks
    return mse, k, beta, iter


def grid_search(
    data: np.ndarray | dict,
    data_type: str,
    n_pca: int | None,
    k_range: np.ndarray,
    beta_range: np.ndarray,
    n_iter: int,
    seed: int = 42,
    n_threads: int = 1,
    progress: bool = True,
) -> pd.DataFrame:
    """
    Perform a grid search of MELD hyperparameters using simulations to test performance
    of parameter combinations.

    Args:
        data (np.ndarray | dict): Precomputed normalised and scaled data matrix,  PCA matrix
        or dictionary of precomputed graphs.
        data_type (str): Type of data matrix ('pca', 'scaled' or 'graph').
        n_pca (int | None): Number of PCA components to use.
        k_range (np.ndarray): Search space for k.
        beta_range (np.ndarray): Search space for beta.
        n_iter (int): Number of simulation iterations for each parameter combination.
        seed (int): Random seed.
        n_threads (int): Number of threads.
        progress (bool): Display progress bars.

    Returns:
        pd.DataFrame: A DataFrame containing the results of the grid search, including
        the mean squared error (mse), the value of k, the value of beta, and the seed.
    """
    # Set default n_pca to 50 if n_pca not specified for normalised and scaled data
    if data_type == "scaled" and n_pca is None:
        warnings.warn(
            "n_pca not specified for normalised and scaled data; using default value of 50"
        )
        n_pca = 50

    with Parallel(n_jobs=int(n_threads), return_as="generator") as parallel:
        if data_type == "graph":
            graphs = data
        else:
            # Precompute k-nearest neighbors graphs (this is the most time-consuming step so perform outside of grid search)
            logger.info("Computing k-nearest neighbors graphs")
            graphs = dict(
                tqdm(
                    parallel(
                        delayed(_compute_graph)(data, n_pca, k, seed) for k in k_range
                    ),
                    desc="Progress",
                    total=len(k_range),
                    disable=(n_threads == 1 or not progress),
                )
            )

        # Initialise Benchmarker for each graph
        logger.info("Initialising Benchmarkers")
        bm = dict(
            tqdm(
                (
                    _initialise_benchmarker(graph, k, n_threads, seed)
                    for k, graph in graphs.items()
                ),
                desc="Progress",
                total=len(k_range),
                disable=(n_threads == 1 or not progress),
            )
        )

        # Run grid search using simulations
        logger.info(
            f"Starting grid search using {len(k_range)} k values, {len(beta_range)} beta values, {n_iter} simulations per combination"
        )
        results = list(
            tqdm(
                parallel(
                    delayed(_run_simulation)(k, bm[k], beta, iter)
                    for (k, beta, iter) in itertools.product(
                        k_range, beta_range, range(n_iter)
                    )
                ),
                desc="Progress",
                total=len(k_range) * len(beta_range) * n_iter,
                disable=(n_threads == 1 or not progress),
            )
        )
        logger.info(
            f"Completed grid search using {len(k_range) * len(beta_range) * n_iter} simulations"
        )

    return pd.DataFrame(results, columns=["mse", "k", "beta", "seed"])


# ==============================
# SCRIPT
# ==============================
if __name__ == "__main__":
    # Run the main function
    _main(opt=docopt.docopt(DOC))
