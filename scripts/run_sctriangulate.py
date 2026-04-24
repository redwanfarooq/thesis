#!/bin/env python


"""
Run scTriangulate.

Requires a MuData or AnnData object containing precomputed cluster labels/reference-based
annotations (in obs columns), UMAP embedding (in obsm slot named 'X_umap'), and optionally,
doublet scores (in obs column named 'doublet_scores').
"""


# ==============================
# MODULES
# ==============================
import anndata as ad
import muon as mu
import sctriangulate as st
import warnings
import docopt
import gc
from loguru import logger


# ==============================
# COMMAND LINE OPTIONS
# ==============================
# Define options
DOC = """
Run scTriangulate.

Requires a MuData or AnnData object containing precomputed cluster labels/reference-based
annotations (in obs columns), UMAP embedding (in obsm slot named 'X_umap'), and optionally,
doublet scores (in obs column named 'doublet_scores').

Usage:
    run_sctriangulate.py --input=<input> --outdir=<outdir> --query=<query> [--counts-layer=<layer>] [options]

Arguments:
    --input=<input>                 Path to input H5AD or H5MU file
    --outdir=<outdir>               Path to the output directory
    --query=<query>                 Semicolon-separated list of query annotations to use for triangulation

Options:
    -t --threads=<threads>          Number of threads to use [default: 1]
    -l --log=<log>                  Path to the log file [default: run_sctriangulate.log]
    --umap-modality=<modality>      Modality to use for UMAP embedding [default: RNA]
    --counts-layer=<layer>          Layer to use for raw counts data
    --species=<species>             Species for the analysis [default: human]
    --criterion=<criterion>         Criterion for selecting artifact genes [default: 2]

    --predict-doublet               Compute doublet scores (will use precomputed doublet scores if 'doublet_scores' column present in obs)
    --compute-metrics-parallel      Compute metrics in parallel
    --compute-shapley-parallel      Compute Shapley values in parallel
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
        "--input",
        "--outdir",
        "--query",
        "--threads",
        "--log",
        "--umap-modality",
        "--counts-layer",
        "--species",
        "--criterion",
        "--predict-doublet",
        "--compute-metrics-parallel",
        "--compute-shapley-parallel",
        "--quiet",
    ]:
        logger.info(f"{key.replace('--', '')}: {opt[key]}")

    # Load input data
    if opt["--input"].endswith(".h5mu"):
        logger.info("Loading MuData object")
        mdata = mu.read_h5mu(opt["--input"])
        adata = ad.concat(
            [mdata[mod].copy() for mod in mdata.mod_names],
            axis=1,
            join="outer",
            merge=None,
            label="modality",
            keys=mdata.mod_names,
        )
        adata.var = adata.var.loc[:, ["modality"]].copy()
        adata.obs = mdata.obs.copy()
        adata.obsm["X_umap"] = mdata[opt["--umap-modality"]].obsm["X_umap"].copy()
        del adata.varm
        del adata.varp
        del adata.raw
        del mdata
        # Force garbage collection to free memory from large MuData object
        gc.collect()
    elif opt["--input"].endswith(".h5ad"):
        logger.info("Loading AnnData object")
        adata = ad.read_h5ad(opt["--input"])
    else:
        raise ValueError("Input file must be in H5AD or H5MU format.")

    # Check for required fields
    assert all([_ in adata.obs.columns for _ in opt["--query"].split(";")]), "Some query annotations not found in obs"
    assert "X_umap" in adata.obsm, "UMAP embedding not found"
    if opt["--counts-layer"] is not None:
        assert opt["--counts-layer"] in adata.layers, f"Counts layer not found"

    # Run scTriangulate
    logger.info("Starting scTriangulate run")
    sctri = st.ScTriangulate(
        dir=opt["--outdir"],
        adata=adata,
        query=opt["--query"].split(";"),
        species=str(opt["--species"]),
        criterion=int(opt["--criterion"]),
        predict_doublet="precomputed" if "doublet_scores" in adata.obs.columns else opt["--predict-doublet"],
        verbose=2,
    )
    sctri.lazy_run(
        layer=opt["--counts-layer"],
        compute_metrics_parallel=opt["--compute-metrics-parallel"],
        compute_shapley_parallel=opt["--compute-shapley-parallel"],
        cores=int(opt["--threads"]),
    )
    logger.success(f"scTriangulate run completed successfully; results saved to {opt['--outdir']}")


# ==============================
# SCRIPT
# ==============================
if __name__ == "__main__":
    # Run the main function
    _main(opt=docopt.docopt(DOC))