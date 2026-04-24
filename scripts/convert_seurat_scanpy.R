#!/bin/env -S Rscript --vanilla


# ==============================
# COMMAND LINE OPTIONS
# ==============================
# Define options
DOC <- "
Convert Seurat object to AnnData/MuData and save to HDF5 format

Usage:
  convert_seurat_scanpy.R --input=<input> --outdir=<outdir> [options]

Arguments:
  REQUIRED
  --input=<input>                                    Path to RDS or QS file containing Seurat object
  --outdir=<outdir>                                  Path to the output directory (file will be named based on input)

  OPTIONAL
  -t --threads=<int>                                 Number of threads [default: 1]
  -l --log=<file>                                    Path to log file [default: convert_seurat_scanpy.log]
  --assays=<assay[;assay...]>                        Assays to include
  --X=<layer>                                        Default layer to include [default: data]
  --layers=<layer[;layer...]>                        Additional layers to include [default: counts;data]
  --obsm=<[assay=]reduction[;[assay=]reduction...]>  Reductions to include (for multiple modalities, specify per modality e.g. RNA=pca;ADT=apca)
  --obsm-global=<reduction[;reduction...]>           Global reductions to include (if multiple modalities present)

Options:
  -h --help                                          Show this screen
  -q --quiet                                         Do not print logging messages to console
"

# Parse options
opt <- docopt::docopt(DOC)

# Logging options
if (!dir.exists(dirname(opt[["log"]]))) dir.create(dirname(opt[["log"]]), recursive = TRUE)
logger::log_appender(logger::appender_file(opt[["log"]]), index = 1)
logger::log_layout(logger::layout_glue, index = 1)
if (!opt[["quiet"]]) {
  logger::log_appender(logger::appender_console, index = 2)
  logger::log_layout(logger::layout_glue_colors, index = 2)
}
logger::log_warnings()
logger::log_errors()

# Extract parameters from command line options
logger::log_info("Parsing command line arguments")
params <- lapply(
  opt,
  function(x) {
    if (is.character(x)) {
      split <- stringr::str_split_1(string = x, pattern = ";")
      if (any(grepl(pattern = "=", x = split))) {
        if (!all(grepl(pattern = "=", x = split))) stop("Some values in 'obsm' are missing names; must provide either all or none.")
        names <- stringr::str_trim(stringr::str_split_i(string = split, pattern = "=", i = 1))
        values <- stringr::str_trim(stringr::str_split_i(string = split, pattern = "=", i = 2))
        names(values) <- names
        return(as.list(values))
      }
      return(split)
    }
    return(x)
  }
)

# Log options
for (x in grep(pattern = "^--|^help", names(params), value = TRUE, invert = TRUE)) {
  y <- if (!is.null(names(params[[x]])) && length(names(params[[x]]))) {
    mapply(
      function(k, v) paste(k, v, sep = "="),
      names(params[[x]]),
      params[[x]],
      SIMPLIFY = TRUE
    )
  } else {
    params[[x]]
  }
  logger::log_info("{x}: {paste(y, collapse = ';')}")
}


# ==============================
# SETUP
# ==============================
logger::log_info("Initialising")

suppressPackageStartupMessages({
  library(reticulate)
  library(Seurat)
  library(BPCells)
})

if (!dir.exists(params$outdir)) dir.create(params$outdir, recursive = TRUE)


# ==============================
# FUNCTIONS
# ==============================
write.h5ad <- function(object, file, assay = NULL, X = "data",
                       layers = c("counts", "data"),
                       obsm = NULL,
                       verbose = TRUE) {
  if (verbose) message("Loading Python libraries")
  ad <- reticulate::import("anndata")

  if (!is(object, "Seurat")) stop("The input object must be a Seurat object")

  if (is.null(assay)) assay <- SeuratObject::DefaultAssay(object)
  if (!assay %in% SeuratObject::Assays(object)) stop("The following assay was not found: ", assay)

  X <- match.arg(X, choices = c("counts", "data"))
  if (!X %in% SeuratObject::Layers(object)) stop("The following layer was not found in assay '", assay, "': ", X)

  if (!is.null(layers)) {
    layers <- match.arg(layers, choices = c("counts", "data"), several.ok = TRUE)
    layers <- layers[layers != X]
    invalid.layers <- !layers %in% SeuratObject::Layers(object)
    layers <- layers[!invalid.layers]
    if (any(invalid.layers) && verbose) warning("The following layers were not found in assay '", assay, "': ", paste(layers[invalid.layers], collapse = ", "))
  }

  if (!is.null(obsm)) {
    if (is.vector(obsm) && !is.list(obsm)) {
      obsm <- as.list(obsm)
    }
    if (is.null(names(obsm))) {
      names(obsm) <- obsm
    }

    invalid.obsm <- !obsm %in% SeuratObject::Reductions(object)
    obsm <- obsm[!invalid.obsm]
    if (any(invalid.obsm) && verbose) warning("The following reductions were not found in the Seurat object: ", paste(obsm[invalid.obsm], collapse = ", "))
  }

  if (verbose) message("Extracting data from Seurat object")
  X <- SeuratObject::GetAssayData(object, assay = assay, layer = X) |> Matrix::t()
  obs <- object[[]]
  var <- object[[assay]][[]]
  if (length(obsm)) {
    obsm <- sapply(obsm, function(x) SeuratObject::Embeddings(object, reduction = x), simplify = FALSE)
  } else {
    obsm <- NULL
  }
  if (length(layers)) layers <- setNames(nm = layers) |> lapply(function(x) SeuratObject::GetAssayData(object, assay = assay, layer = x) |> Matrix::t())

  if (verbose) message("Creating AnnData object")
  adata <- ad$AnnData(
    X = X,
    obs = obs,
    var = var,
    obsm = obsm,
    layers = layers
  )

  if (verbose) message("Writing AnnData object to file")
  adata$write(file, compression = "gzip")

  invisible(file)
}


write.h5mu <- function(object, file,
                       assays = NULL,
                       X = "data",
                       layers = c("counts", "data"),
                       obsm = NULL,
                       obsm.global = NULL,
                       verbose = TRUE) {
  if (verbose) message("Loading Python libraries")
  ad <- reticulate::import("anndata")
  mu <- reticulate::import("mudata")
  mu$set_options(pull_on_update = FALSE)

  if (!is(object, "Seurat")) stop("The input object must be a Seurat object")

  available.assays <- SeuratObject::Assays(object)
  if (is.null(assays)) {
    assays <- available.assays
  } else {
    missing.assays <- !assays %in% available.assays
    if (any(missing.assays)) {
      stop("The following assays were not found: ", paste(assays[missing.assays], collapse = ", "))
    }
  }

  if (length(X) == 1 && is.null(names(X))) {
    X <- setNames(rep(X, length(assays)), assays)
  }

  if (is.character(layers) || (is.list(layers) && is.null(names(layers)))) {
    layers <- setNames(rep(list(layers), length(assays)), assays)
  }

  if (is.character(obsm) || (is.list(obsm) && is.null(names(obsm)))) {
    obsm <- setNames(rep(list(obsm), length(assays)), assays)
  }

  adata.list <- list()

  if (verbose) message("Extracting data from Seurat object")
  for (assay in assays) {
    if (verbose) message("Processing assay: ", assay)

    X.layer <- if (assay %in% names(X)) X[[assay]] else "data"
    X.layer <- match.arg(X.layer, choices = c("counts", "data"))

    if (!X.layer %in% SeuratObject::Layers(object, assay = assay)) {
      stop("The following layer was not found in assay '", assay, "': ", X.layer)
    }

    assay.layers <- NULL
    if (assay %in% names(layers) && !is.null(layers[[assay]])) {
      assay.layers <- match.arg(layers[[assay]], choices = c("counts", "data"), several.ok = TRUE)
      assay.layers <- assay.layers[assay.layers != X.layer]
      invalid.layers <- !assay.layers %in% SeuratObject::Layers(object, assay = assay)
      assay.layers <- assay.layers[!invalid.layers]
      if (any(invalid.layers) && verbose) {
        warning(
          "The following layers were not found in assay '", assay, "': ",
          paste(assay.layers[invalid.layers], collapse = ", ")
        )
      }
    }

    X.data <- SeuratObject::GetAssayData(object, assay = assay, layer = X.layer) |> Matrix::t()

    layers.data <- NULL
    if (length(assay.layers)) {
      layers.data <- setNames(nm = assay.layers) |>
        lapply(function(x) SeuratObject::GetAssayData(object, assay = assay, layer = x) |> Matrix::t())
    }

    var <- object[[assay]][[]]

    obsm.data <- NULL
    if (!is.null(obsm) && assay %in% names(obsm)) {
      assay.obsm <- obsm[[assay]]

      if (is.vector(assay.obsm) && !is.list(assay.obsm)) {
        assay.obsm <- as.list(assay.obsm)
      }
      if (is.null(names(assay.obsm))) {
        names(assay.obsm) <- assay.obsm
      }

      if (length(assay.obsm)) {
        invalid.obsm <- !assay.obsm %in% SeuratObject::Reductions(object)
        assay.obsm <- assay.obsm[!invalid.obsm]
        if (any(invalid.obsm) && verbose) {
          warning(
            "The following reductions were not found for assay '", assay, "': ",
            paste(assay.obsm[invalid.obsm], collapse = ", ")
          )
        }
        if (length(assay.obsm)) {
          obsm.data <- sapply(assay.obsm, function(x) SeuratObject::Embeddings(object, reduction = x), simplify = FALSE)
        }
      }
    }

    adata.list[[assay]] <- ad$AnnData(
      X = X.data,
      obs = data.frame(row.names = colnames(object[[assay]])),
      var = var,
      obsm = obsm.data,
      layers = layers.data
    )
  }

  if (verbose) message("Creating MuData object")

  obsm.global.data <- NULL
  if (!is.null(obsm.global)) {
    if (is.vector(obsm.global) && !is.list(obsm.global)) {
      obsm.global <- as.list(obsm.global)
    }
    if (is.null(names(obsm.global))) {
      names(obsm.global) <- obsm.global
    }

    invalid.obsm.global <- !obsm.global %in% SeuratObject::Reductions(object)
    obsm.global <- obsm.global[!invalid.obsm.global]
    if (any(invalid.obsm.global) && verbose) {
      warning(
        "The following global reductions were not found in the Seurat object: ",
        paste(obsm.global[invalid.obsm.global], collapse = ", ")
      )
    }

    if (length(obsm.global)) {
      obsm.global.data <- sapply(obsm.global, function(x) SeuratObject::Embeddings(object, reduction = x), simplify = FALSE)
    }
  }

  mdata <- mu$MuData(adata.list)
  # Assign global obs data directly (identical across all modalities)
  mdata$obs <- object[[]]
  mdata$obsm <- obsm.global.data

  if (verbose) message("Writing MuData object to file")
  mdata$write(file, compression = "gzip")

  invisible(file)
}

# ==============================
# SCRIPT
# ==============================
logger::log_info("Loading input Seurat object")
# Load input
if (grepl(pattern = "\\.qs$", x = params$input, ignore.case = TRUE)) {
  seu <- qs::qread(file = params$input, nthreads = as.integer(params$threads))
} else {
  seu <- readRDS(file = params$input)
}
nm <- tools::file_path_sans_ext(basename(params$input))

# Convert to AnnData/MuData and save to HDF5 format
if (length(SeuratObject::Assays(seu)) == 1 || length(params$assays) == 1) {
  logger::log_info("Converting to AnnData")
  params$obsm <- c(params$obsm, params$obsm_global) |> unique()
  out <- write.h5ad(
    object = seu,
    file = file.path(params$outdir, paste0(nm, ".h5ad")),
    assay = params$assays,
    X = params$X,
    layers = params$layers,
    obsm = params$obsm,
    verbose = !params$quiet
  )
} else if (length(SeuratObject::Assays(seu)) > 1) {
  logger::log_info("Converting to MuData")
  out <- write.h5mu(
    object = seu,
    file = file.path(params$outdir, paste0(nm, ".h5mu")),
    assays = params$assays,
    X = params$X,
    layers = params$layers,
    obsm = params$obsm,
    obsm.global = params$obsm_global,
    verbose = !params$quiet
  )
} else {
  stop("The input Seurat object does not contain any assays. Please check the input file.")
}
logger::log_success("Output path: {out}")
