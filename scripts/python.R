require(reticulate)


#' Convert Seurat object to AnnData and write to h5ad file
#'
#' This function converts a Seurat object to an AnnData object and writes it to an h5ad file.
#'
#' @param object A Seurat object.
#' @param file A character string specifying the path to the output h5ad file.
#' @param assay A character string specifying the assay to use for the AnnData object. Default is `"RNA"`.
#' @param X A character string specifying the layer to use for the AnnData object. Default is `"data"`.
#' @param layers A character vector specifying the additional layers to include in the AnnData object. Default is
#' `c("counts", "data")`.
#' @param obsm A (optionally named) list specifying the embeddings to include in the AnnData object. Default is
#' NULL.
#' @param verbose A logical value indicating whether to print progress messages. Default is `TRUE`.
#'
#' @return The path to the output h5ad file.
#'
#' @importFrom reticulate import
#' @importFrom SeuratObject DefaultAssay Assays GetAssayData Layers Reductions Embeddings
#'
#' @export
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


#' Convert multimodal Seurat object to MuData and write to h5mu file
#'
#' This function converts a multimodal Seurat object to a MuData object and writes it to an h5mu file.
#' Each assay in the Seurat object becomes a separate modality in the MuData object. The global
#' cell metadata from the Seurat object is stored in the MuData `.obs` slot (not in individual
#' AnnData objects), while feature metadata is stored in the `.var` slot of each AnnData object.
#' Assay-specific reductions are stored in the `.obsm` slot of individual AnnData objects, while
#' global reductions are stored in the MuData `.obsm` slot.
#'
#' @param object A Seurat object with multiple assays.
#' @param file A character string specifying the path to the output h5mu file.
#' @param assays A character vector specifying the assays to include. If NULL, all assays will be included.
#' @param X A named character vector specifying the layer to use for each assay's X matrix, or a single character string
#' to use the same layer for all assays. Names should match assay names when using named vector.
#' Example: `c("RNA" = "data", "ATAC" = "data")` or `"data"`.
#' @param layers A named list specifying additional layers to include for each assay, or a character vector to use
#' the same layers for all assays. Example: `list("RNA" = c("counts", "data"), "ATAC" = c("counts", "data"))` or `c("counts", "data")`.
#' @param obsm A named list specifying the assay-specific embeddings to include in each AnnData object's obsm slot.
#' Should be structured as `list("assay.name" = c("embedding.name" = "reduction.name"))` or
#' `list("assay.name" = list("embedding.name" = "reduction.name"))`. Default is `NULL`.
#' @param obsm.global A named list specifying the multimodal embeddings to include in the MuData object's obsm slot.
#' These are shared across all modalities. Example: `list(pca = "pca", umap = "umap")`.
#' @param verbose A logical value indicating whether to print progress messages. Default is `TRUE`.
#'
#' @return The path to the output h5mu file.
#'
#' @importFrom reticulate import
#' @importFrom SeuratObject Assays GetAssayData Layers Reductions Embeddings
#'
#' @export
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
  mdata$obs <- object[[]]
  mdata$obsm <- obsm.global.data

  if (verbose) message("Writing MuData object to file")
  mdata$write(file, compression = "gzip")

  invisible(file)
}


#' Convert AnnData (h5ad) file to Seurat object
#'
#' This function reads an AnnData h5ad file and converts it to a Seurat object.
#' Supports dense and sparse matrices, obs/var metadata, and selected dimension reductions from obsm.
#'
#' @param file Path to the h5ad file.
#' @param x Character string specifying which layer to use for the main data matrix (default: "data").
#' @param raw Character string specifying which layer to use for the raw counts matrix (default: "counts").
#' @param assay.name Name for the Seurat assay (default: "RNA").
#' @param obsm Named list specifying which dimension reductions to import from obsm (e.g. list(pca = "X_pca")).
#'
#' @return A Seurat object containing the data and metadata from the h5ad file.
#'
#' @importFrom hdf5r H5File
#' @importFrom Matrix sparseMatrix
#' @importFrom Seurat CreateSeuratObject SetAssayData CreateDimReducObject
#'
#' @export
h5ad.to.seurat <- function(
    file,
    x = "data",
    raw = "counts",
    assay.name = "RNA",
    obsm = NULL) {
  # Input checks
  if (!is.character(file) || length(file) != 1 || !file.exists(file)) {
    stop("File must be a valid path to an existing h5ad file.")
  }
  if (!is.character(assay.name) || length(assay.name) != 1) {
    stop("'assay.name' must be a single string.")
  }
  valid.layers <- c("counts", "data", "scale.data")
  if (!is.null(x) && !(x %in% valid.layers)) {
    stop("'x' must be one of: ", paste(valid.layers, collapse = ", "))
  }
  if (!is.null(raw) && !(raw %in% valid.layers)) {
    stop("'raw' must be one of: ", paste(valid.layers, collapse = ", "))
  }

  # Helper: load matrix (dense or sparse)
  load.h5ad.matrix <- function(h5.obj, dimnames = NULL) {
    if (inherits(h5.obj, "H5D")) {
      mat <- h5.obj[, ]
      mat <- t(mat)
    } else if (inherits(h5.obj, "H5Group")) {
      group.names <- names(h5.obj)
      if (all(c("indptr", "indices", "data") %in% group.names)) {
        indptr <- h5.obj[["indptr"]][]
        indices <- h5.obj[["indices"]][]
        data <- h5.obj[["data"]][]
        shape <- h5.obj$attr_open("shape")$read()
        mat <- Matrix::sparseMatrix(
          i = indices + 1,
          p = indptr,
          x = data,
          dims = shape,
          dimnames = dimnames,
          index1 = TRUE
        )
      } else {
        stop("Sparse matrix group missing required datasets or shape attribute.")
      }
    } else {
      stop("Input is not a valid HDF5 dataset or group.")
    }
    if (!is.null(dimnames)) {
      dimnames(mat) <- dimnames
    }
    return(mat)
  }

  # Helper: load table from AnnData obs/var, including categorical columns
  load.h5ad.table <- function(h5.group) {
    col.names <- names(h5.group)
    columns <- list()
    for (n in col.names) {
      obj <- h5.group[[n]]
      if (inherits(obj, "H5D")) {
        columns[[n]] <- obj[]
      } else if (inherits(obj, "H5Group")) {
        if (all(c("codes", "categories") %in% names(obj))) {
          codes <- obj[["codes"]][]
          categories <- obj[["categories"]][]
          col <- rep(NA_character_, length(codes))
          valid <- codes != -1
          col[valid] <- categories[codes[valid] + 1]
          columns[[n]] <- col
        } else {
          warning(sprintf("Group '%s' does not have both 'codes' and 'categories'. Skipping.", n))
        }
      } else {
        warning(sprintf("Column '%s' is not H5D or H5Group. Skipping.", n))
      }
    }
    df <- as.data.frame(columns, stringsAsFactors = FALSE)
    if ("X_index" %in% names(df)) {
      rownames(df) <- df[["X_index"]]
      df[["X_index"]] <- NULL
    }
    return(df)
  }

  # Open h5ad file
  h5.file <- hdf5r::H5File$new(filename = file, mode = "r")
  on.exit(h5.file$close_all())

  has.group <- function(name) name %in% names(h5.file)
  is.dense <- function(obj) inherits(obj, "H5D")

  # Load obs and var tables
  if (!has.group("obs") || !has.group("var")) {
    stop("h5ad file must contain 'obs' and 'var' groups.")
  }
  obs <- load.h5ad.table(h5.file[["obs"]])
  var <- load.h5ad.table(h5.file[["var"]])
  raw.exists <- has.group("raw")
  var.raw <- if (raw.exists) load.h5ad.table(h5.file[["raw/var"]]) else NULL

  # Load X and raw/X matrices
  x.exists <- has.group("X")
  x.obj <- if (x.exists) h5.file[["X"]] else NULL
  raw.obj <- if (raw.exists) h5.file[["raw/X"]] else NULL

  if (is.null(x.obj)) stop("h5ad file missing 'X' dataset.")

  x.mat <- t(load.h5ad.matrix(x.obj, dimnames = list(rownames(obs), rownames(var))))
  raw.mat <- if (!is.null(raw.obj) && !is.null(var.raw)) {
    t(load.h5ad.matrix(raw.obj, dimnames = list(rownames(obs), rownames(var.raw))))
  } else {
    NULL
  }

  # Assign matrices to layers
  mat <- list(counts = NULL, data = NULL, scale.data = NULL)
  if (!is.null(x)) mat[[x]] <- x.mat
  if (!is.null(raw) && !is.null(raw.mat)) mat[[raw]] <- raw.mat

  # Fill unset layers with defaults only if x and raw are NULL
  if (is.null(x) && is.null(raw)) {
    if (raw.exists && !is.null(raw.mat)) {
      if (is.dense(x.obj)) {
        mat$data <- raw.mat
        mat$scale.data <- x.mat
      } else {
        mat$counts <- raw.mat
        mat$data <- x.mat
      }
    } else {
      mat$data <- x.mat
    }
  }

  # Create Seurat object
  suppressWarnings({
    seurat.obj <- Seurat::CreateSeuratObject(
      counts = mat$counts,
      data = mat$data,
      meta.data = obs,
      assay = assay.name
    )
    if (!is.null(mat$scale.data)) {
      seurat.obj <- Seurat::SetAssayData(
        object = seurat.obj,
        assay = assay.name,
        layer = "scale.data",
        new.data = mat$scale.data
      )
    }
    if (!is.null(var.raw)) {
      combined.var <- var
      missing.rows <- setdiff(rownames(var.raw), rownames(var))
      if (length(missing.rows) > 0) {
        na.df <- as.data.frame(matrix(NA, nrow = length(missing.rows), ncol = ncol(var)))
        colnames(na.df) <- colnames(var)
        rownames(na.df) <- missing.rows
        combined.var <- rbind(combined.var, na.df)
      }
      for (col in intersect(colnames(var.raw), colnames(var))) {
        combined.var[rownames(var.raw), col] <- var.raw[, col]
      }
      new.cols <- setdiff(colnames(var.raw), colnames(var))
      if (length(new.cols) > 0) {
        combined.var[rownames(var.raw), new.cols] <- var.raw[, new.cols]
      }
      combined.var <- combined.var[rownames(seurat.obj[[assay.name]]), , drop = FALSE]
      seurat.obj[[assay.name]][[]] <- combined.var
    } else {
      seurat.obj[[assay.name]][[]] <- var
    }
  })

  # Dimension reductions
  if (!is.null(obsm) && has.group("obsm")) {
    obsm.group <- h5.file[["obsm"]]
    for (dr.name in names(obsm)) {
      obsm.key <- obsm[[dr.name]]
      if (obsm.key %in% names(obsm.group)) {
        dr.mat <- t(obsm.group[[obsm.key]][, ])
        if (nrow(dr.mat) != nrow(obs)) {
          warning(sprintf("Dimension reduction '%s' has %d rows, expected %d (number of cells). Skipping.", dr.name, nrow(dr.mat), nrow(obs)))
          next
        }
        colnames(dr.mat) <- paste0(dr.name, "_", seq_len(ncol(dr.mat)))
        rownames(dr.mat) <- rownames(obs)
        suppressWarnings({
          seurat.obj[[dr.name]] <- Seurat::CreateDimReducObject(
            embeddings = dr.mat,
            key = paste0(dr.name, "_"),
            assay = assay.name
          )
        })
      } else {
        warning(sprintf("obsm key '%s' not found for dimension reduction '%s'.", obsm.key, dr.name))
      }
    }
  }

  return(seurat.obj)
}


#' Run MELD analysis
#'
#' This function is an R wrapper for Python functions to run the MELD algorithm to compute condition-specific relative
#' likelihoods for each cell per replicate.
#'
#' @param data A matrix of single-cell data (either normalised counts, low-dimensional embeddings or precomputed graph
#' adjacency/affinity matrix). Rows correspond to cells and columns correspond to features.
#' @param cond A factor indicating the condition for each row in `data`.
#' @param rep An optional vector specifying the replicate for each row in `data`. If not provided, all rows will
#' be considered as originating from a single replicate.
#' @param n_pca An integer specifying the number of PCA components to use for calculating neighborhoods. Default is NULL,
#' which will use the original data.
#' @param knn An integer specifying the number of nearest neighbors to consider when constructing the kNN graph.
#' @param beta An integer specifying the beta parameter for the MELD algorithm.
#' @param simplify A logical value indicating whether to simplify the output by keeping only the likelihoods for
#' the first level of `cond`.
#'
#' @return A data frame containing condition-specific relative likelihoods for each cell per replicate.
#'
#' @details This function performs the following steps:
#'   1. Computes the kNN graph using the input data and the specified value of k.
#'   2. Fits the MELD model to the graph using the sample labels.
#'   3. Normalizes the sample densities obtained from the MELD model.
#'   4. Optionally simplifies the output by keeping only the likelihoods for the first level of the condition factor.
#'   5. Returns the sample likelihoods as a data frame.
#'
#' @references
#' Burkhardt, D. B., Stanley, J. S, Tong, A. Perdigoto, A. L., Gigante, S., Herold, K. C., Wolf, G., Giraldez, A. J.,
#' van Dijk, D., and Krishnaswamy, S. (2021). Quantifying the effect of experimental perturbations at single-cell
#' resolution. Nature Biotechnology, 39, 619-629.
#'
#' @importFrom reticulate import
#'
#' @export
run.meld <- function(data,
                     cond,
                     rep = NULL,
                     n_pca = NULL,
                     knn = 5,
                     beta = 60,
                     simplify = FALSE,
                     verbose = TRUE) {
  if (!is.null(n_pca)) if (n_pca > ncol(data)) stop("Number of PCA components exceeds the number of columns in the input matrix: ", ncol(data))
  if (!is.factor(cond)) stop("'cond' must be a factor")
  knn <- as.integer(knn)
  beta <- as.integer(beta)

  if (verbose) message("Loading Python libraries")
  pd <- reticulate::import("pandas")
  gt <- reticulate::import("graphtools")
  meld <- reticulate::import("meld")

  if (nrow(data) == ncol(data)) {
    if (verbose) message("Using precomputed graph")
    graph <- gt$Graph(data, precomputed = "adjacency", use_pygsp = TRUE)
  } else {
    if (verbose) message(glue::glue("Computing kNN graph with k = {knn}"))
    graph <- gt$Graph(data, n_pca = if (!is.null(n_pca)) as.integer(n_pca) else n_pca, knn = knn, use_pygsp = TRUE)
  }

  if (is.null(rep)) rep <- "R1"
  labels <- pd$Series(paste(rep, cond, sep = "_"))

  if (verbose) message(glue::glue("Running MELD with beta = {beta}"))
  meld.op <- meld$MELD(beta = beta)
  sample.densities <- meld.op$fit_transform(graph, sample_labels = labels)

  if (verbose) message("Normalizing sample densities")
  replicates <- unique(rep)
  groups <- lapply(replicates, function(x) grep(pattern = paste0("^", x, "_"), x = colnames(sample.densities), value = TRUE))
  sample.likelihoods <- lapply(groups, function(cols) apply(sample.densities[, cols], MARGIN = 1, function(x) x / sum(x)))
  sample.likelihoods <- do.call(rbind, sample.likelihoods) |>
    t() |>
    as.data.frame()
  rownames(sample.likelihoods) <- rownames(data)

  if (simplify) {
    sel <- levels(cond)[1]
    if (verbose) message(glue::glue("Getting likelihoods for condition '{sel}'"))
    sample.likelihoods <- sample.likelihoods[, grep(pattern = paste0("_", sel, "$"), x = colnames(sample.likelihoods))]
    colnames(sample.likelihoods) <- gsub(pattern = paste0("_", sel, "$"), replacement = "", x = colnames(sample.likelihoods))
  }

  return(sample.likelihoods)
}
