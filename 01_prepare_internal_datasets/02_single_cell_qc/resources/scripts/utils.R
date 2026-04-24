# ==============================
# Helper functions
# ==============================


require(Rcpp)


# Explicitly define pipe and set if null operator
`%>%` <- magrittr::`%>%`
`%||%` <- rlang::`%||%`


# C++ function for efficiently reading BarCounter CSV file to sparse matrix
Rcpp::sourceCpp("read_barcounter_csv.cpp")


#' Modified implementation of CellRanger OrdMag algorithm
#'
#' Finds the `q`th quantile of `x` and returns a threshold `r`-fold less
#' than this value.
#'
#' @param x Numeric vector of log10-transformed counts.
#' @param q Numeric scalar between 0 and 1 (default 0.99).
#' @param r Numeric scalar greater than 1 (default 10).
#'
#' @returns Numeric scalar.
.ordmag <- function(x, q = 0.99, r = 10) {
  return(quantile(x, q, na.rm = TRUE) - log10(r))
}


#' Update R Markdown params list
#'
#' Recursively searches nested named list of default parameters for matching entries
#' in a user-defined list. If a match is found, the default value is replaced with
#' the user-defined value.
#'
#' @param user Named list containing user-defined parameters.
#' @param default Named list containing default parameters (named `params` by default
#'  in R Markdown documents).
#'
#' @returns A nested named list with the same structure as `default`.
#'
#' @export
update.params <- function(user, default = params) {
  if (!all(is.list(user), !is.null(names(user)))) stop("'user' must be a named list")
  if (!all(is.list(default), !is.null(names(default)))) stop("'default' must be a named list")

  new <- vector("list", length(default)) %>% setNames(names(default))
  for (x in names(default)) {
    if (x %in% names(user)) {
      if (is.list(default[[x]])) {
        new[[x]] <- update.params(user[[x]], default[[x]])
      } else {
        new[[x]] <- user[[x]]
      }
    } else {
      new[[x]] <- default[[x]]
    }
  }

  return(new)
}


#' Prepend index to a string
#'
#' Utility function for generating dynamic numbered lists in Markdown. Searches for
#' a variable named `list.index` in the calling environment which is used to keep a
#' running count of the index (automatically incremented with each call).
#'
#' @param x Character scalar specifying list item.
#'
#' @returns A character scalar of `x` prepended with an index.
#'
#' @export
prepend.index <- function(x) {
  i <- get("list.index", envir = parent.frame()) %>% as.integer()
  str <- glue::glue("{i}. {x}")
  list.index <<- i + 1
  return(str)
}


#' Get 10x count matrix
#'
#' Loads 10x count matrix (in HDF5 format).
#'
#' @param file Path to file.
#' @param version Character scalar. 10x HDF5 version ('auto', 'v2', or 'v3').
#' @param cells Character vector. If specified, will subset to matching cell
#' barcodes.
#' @param type Character vector. If specified, will subset to features with matching
#' values in feature type field (v3 only).
#' @param remove.suffix Logical scalar (default `TRUE`). Remove '-1' suffix automatically
#' appended to cell barcodes by Cell Ranger.
#' @param group Optional character scalar. Group name in HDF5 file containing count matrix.
#' If not specified, will use 'matrix' for v3 and the only group for v2.
#'
#' @returns A sparse matrix of counts with features as row names and cell barcodes
#' as column names. If feature type not specified returns a named list of sparse matrices
#' (one per feature type detected).
#'
#' @importFrom hdf5r H5File existsGroup
#'
#' @export
get.10x.h5 <- function(file,
                       version = c("auto", "v2", "v3"),
                       cells = NULL,
                       type = NULL,
                       remove.suffix = TRUE,
                       group = NULL) {
  if (!file.exists(file)) stop(file, " does not exist")

  infile <- hdf5r::H5File$new(filename = file, mode = "r")

  version <- match.arg(version)
  if (version == "auto") {
    version <- if (hdf5r::existsGroup(infile, "matrix")) "v3" else "v2"
    message("Detected 10x HDF5 version: ", version)
  }
  features <- if (version == "v2") "genes" else "features/id"

  if (is.null(group)) {
    if (version == "v3") {
      group <- "matrix"
    } else {
      group <- names(infile)
      if (length(group) > 1) stop("Multiple groups detected: ", paste(group, collapse = ", "), ". Please specify 'group' parameter.")
    }
  }
  counts <- infile[[paste(group, "data", sep = "/")]]
  indices <- infile[[paste(group, "indices", sep = "/")]]
  indptr <- infile[[paste(group, "indptr", sep = "/")]]
  shape <- infile[[paste(group, "shape", sep = "/")]]
  features <- infile[[paste(group, features, sep = "/")]]
  barcodes <- infile[[paste(group, "barcodes", sep = "/")]]
  matrix <- Matrix::sparseMatrix(
    i = indices[],
    p = indptr[],
    x = as.numeric(counts[]),
    dims = shape[],
    dimnames = list(features[], barcodes[]),
    index1 = FALSE,
    repr = "C"
  )

  if (remove.suffix) colnames(matrix) <- gsub(pattern = "-1$", replacement = "", colnames(matrix))

  if (!is.null(cells)) matrix <- matrix[, match(cells, colnames(matrix))]
  if (any(is.na(colnames(matrix)))) {
    warning(sum(is.na(colnames(matrix))), " barcode(s) in 'cells' not present in count matrix")
    matrix <- matrix[, !is.na(colnames(matrix))]
  }

  if (version == "v3") {
    feature.type <- infile[["matrix/features/feature_type"]]
    if (!is.null(type)) {
      matrix <- matrix[feature.type[] %in% type, ]
    } else {
      if (length(unique(feature.type[])) > 1) message("Multiple feature types detected: ", paste(unique(feature.type[]), collapse = ", "), ".")
      message("Returning list of matrices; please specify 'type' parameter to return a single matrix.")
      matrix <- lapply(
        unique(feature.type[]),
        function(type, matrix, feature.type) matrix[grep(pattern = type, x = feature.type), ],
        matrix = matrix,
        feature.type = feature.type[]
      ) |>
        setNames(unique(feature.type[]))
    }

    return(matrix)
  }
}


#' Get 10x count matrix
#'
#' Loads 10x count matrix (in matrix market format).
#'
#' @param file Path to file.
#' @param cells Character vector. If specified, will subset to matching cell
#'  barcodes.
#' @param type Character vector. If specified, will subset to features with matching
#'  values in column 3 of 10x features TSV file.
#' @param remove.suffix Logical scalar (default `TRUE`). Remove '-1' suffix automatically
#'  appended to cell barcodes by Cell Ranger.
#'
#' @returns A sparse matrix of counts with features as row names and cell barcodes
#'  as column names.
#'
#' @importFrom Matrix readMM
#'
#' @export
get.10x.matrix <- function(file,
                           cells = NULL,
                           type = NULL,
                           remove.suffix = TRUE) {
  if (!file.exists(file)) stop(file, " does not exist")

  matrix <- Matrix::readMM(file) %>% as("CsparseMatrix")
  barcodes <- get.10x.barcodes(file.path(dirname(file), "barcodes.tsv.gz"), remove.suffix = remove.suffix)
  features <- get.10x.features(file.path(dirname(file), "features.tsv.gz"))
  colnames(matrix) <- rownames(barcodes)
  rownames(matrix) <- rownames(features)
  if (!is.null(cells)) matrix <- matrix[, match(cells, colnames(matrix))]
  if (any(is.na(colnames(matrix)))) {
    warning(sum(is.na(colnames(matrix))), " barcode(s) in 'cells' not present in count matrix")
    matrix <- matrix[, !is.na(colnames(matrix))]
  }
  if (!is.null(type)) matrix <- matrix[features$Type %in% type, ]

  return(matrix)
}


#' Get 10x cell barcodes
#'
#' Loads 10x cell barcodes from TSV.
#'
#' @inheritParams get.10x.matrix
#'
#' @returns A data frame with the following columns:
#'  * Barcode: cell barcode
#'
#' @importFrom vroom vroom
#'
#' @export
get.10x.barcodes <- function(file, remove.suffix = FALSE) {
  if (!file.exists(file)) stop(file, " does not exist")

  df <- vroom::vroom(
    file,
    delim = "\t",
    col_names = "Barcode",
    show_col_types = FALSE
  ) %>%
    as.data.frame()
  if (remove.suffix) df$Barcode <- gsub(pattern = "-1$", replacement = "", df$Barcode)
  rownames(df) <- df$Barcode

  return(df)
}


#' Get 10x feature metadata
#'
#' Loads 10x feature metadata from TSV.
#'
#' @inheritParams get.10x.matrix
#'
#' @returns A data frame with the following columns:
#'  * ID: feature ID
#'  * Symbol: feature symbol
#'  * Type: feature type
#'  * Chr: chromosome
#'  * Start: genomic start position
#'  * End: genomic end position
#'
#' @importFrom vroom vroom
#'
#' @export
get.10x.features <- function(file, type = NULL) {
  if (!file.exists(file)) stop(file, " does not exist")

  df <- vroom::vroom(
    file,
    delim = "\t",
    col_names = c("ID", "Symbol", "Type", "Chr", "Start", "End"),
    show_col_types = FALSE
  ) %>%
    as.data.frame()
  rownames(df) <- df$ID
  if (!is.null(type)) df <- df[df$Type %in% type, ]

  return(df)
}


#' Get BarCounter count matrix
#'
#' Loads BarCounter count matrix (in CSV format).
#'
#' @param file Path to file.
#' @param cells Character vector. If specified, will subset to matching cell
#'  barcodes (after adding suffix if `add.suffix = TRUE`).
#' @param features.pattern Regular expression string. If specified, with subset
#'  to features with matching names.
#' @param include.total Logical scalar (default `FALSE`). Include column containing
#'  total UMI count per barcode.
#' @param add.suffix Logical scalar (default `FALSE`). Add '-1' suffix to match
#'  cell barcodes from Cell Ranger.
#'
#' @returns A sparse matrix of counts with features as row names and cell barcodes
#'  as column names.
#'
#' @export
get.barcounter.matrix <- function(file,
                                  cells = NULL,
                                  features.pattern = NULL,
                                  include.total = FALSE,
                                  add.suffix = FALSE) {
  if (!file.exists(file)) stop(file, " does not exist")

  # Efficiently read BarCounter CSV file into sparse matrix without loading dense matrix into memory
  sparse <- read_barcounter_csv(file)
  matrix <- Matrix::sparseMatrix(
    i = sparse$i,
    p = sparse$p,
    x = sparse$x,
    dims = c(length(sparse$features), length(sparse$barcodes)),
    dimnames = list(sparse$features, sparse$barcodes),
    index1 = FALSE,
    repr = "C"
  )
  if (add.suffix) colnames(matrix) <- paste(colnames(matrix), "1", sep = "-")
  if (!include.total) matrix <- matrix[rownames(matrix) != "total", ]
  if (!is.null(cells)) matrix <- matrix[, match(cells, colnames(matrix))]
  if (any(is.na(colnames(matrix)))) {
    warning(sum(is.na(colnames(matrix))), " barcode(s) in 'cells' not present in count matrix")
    matrix <- matrix[, !is.na(colnames(matrix))]
  }
  if (!is.null(features.pattern)) matrix <- matrix[grepl(pattern = features.pattern, x = rownames(matrix)), ]

  return(matrix)
}


#' Add cell metadata to SingleCellExperiment or Seurat object
#'
#' Merge or replace cell metadata contained in a SingleCellExperiment or
#' Seurat object with values in a data frame or data frame-like object.
#'
#' @param x A SingleCellExperiment or Seurat object.
#' @param metadata A data frame or data frame-like object.
#' @param replace Character scalar:
#'  * 'matching' (default): Replace metadata columns in `x` with
#'    matching columns in `metadata`.
#'  * 'all': Replace all metadata columns in `x` with columns in
#'    `metadata`.
#'  * 'none': Do not replace any metadata columns in `x`.
#'
#' @returns An object of the same class as `x` with updated cell metadata.
#'
#' @importFrom SummarizedExperiment colData
#'
#' @export
add.cell.metadata <- function(x,
                              metadata,
                              replace = c("matching", "all", "none")) {
  if (!any(is(x, "SingleCellExperiment"), is(x, "Seurat"))) {
    stop("x must be an object of class 'SingleCellExperiment' or 'Seurat'")
  }
  replace <- match.arg(replace)

  if (is(x, "SingleCellExperiment")) {
    current <- SummarizedExperiment::colData(x)
  } else if (is(x, "Seurat")) {
    current <- x@meta.data
  }

  if (!setequal(rownames(current), rownames(metadata))) stop("x and metadata must contain the same cell barcodes")
  metadata <- metadata[rownames(current), , drop = FALSE] # ensure matching row order

  if (replace == "matching") {
    current <- current[, !colnames(current) %in% intersect(colnames(current), colnames(metadata))]
  } else if (replace == "all") {
    current <- data.frame(row.names = rownames(current))
  }

  new <- cbind(current, metadata)

  if ((is(x, "SingleCellExperiment"))) {
    SummarizedExperiment::colData(x) <- new %>% as("DataFrame")
  } else if (is(x, "Seurat")) {
    x@meta.data <- new %>% as("data.frame")
  }

  return(x)
}


#' Add feature metadata to SingleCellExperiment or Seurat object
#'
#' Merge or replace feature metadata contained in a SingleCellExperiment or
#' Seurat object with values in a data frame or data frame-like object.
#'
#' @inheritParams add.cell.metadata
#' @param assay Name of assay containing relevant features (if `x` is a Seurat
#'  object). If not provided, will use `DefaultAssay(x)`.
#'
#' @returns An object of the same class as `x` with updated feature metadata.
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom SeuratObject DefaultAssay
#'
#' @export
add.feature.metadata <- function(x,
                                 metadata,
                                 assay = NULL,
                                 replace = c("matching", "all", "none")) {
  if (!any(is(x, "SingleCellExperiment"), is(x, "Seurat"))) {
    stop("x must be an object of class 'SingleCellExperiment' or 'Seurat'")
  }
  replace <- match.arg(replace)

  if (is(x, "SingleCellExperiment")) {
    current <- SummarizedExperiment::rowData(x)
  } else if (is(x, "Seurat")) {
    assay <- ifelse(is.null(assay), SeuratObject::DefaultAssay(x), assay)
    current <- x[[assay]]@meta.features
  }

  if (!setequal(rownames(current), rownames(metadata))) stop("x and metadata must contain the same features")
  metadata <- metadata[rownames(current), , drop = FALSE] # ensure matching row order

  if (replace == "matching") {
    current <- current[, !colnames(current) %in% intersect(colnames(current), colnames(metadata))]
  } else if (replace == "all") {
    current <- data.frame(row.names = rownames(current))
  }

  new <- cbind(current, metadata)

  if ((is(x, "SingleCellExperiment"))) {
    SummarizedExperiment::rowData(x) <- new %>% as("DataFrame")
  } else if (is(x, "Seurat")) {
    x[[assay]]@meta.features <- new %>% as("data.frame")
  }

  return(x)
}


#' Modified implementation of CellRanger empirical expected cell number estimation
#' algorithm
#'
#' Estimates the expected cell number *x* from empirical total UMI counts per barcode
#' by minimising the loss function (OrdMag(*x*) - *x*)^2 / *x* over the range of *x*
#' between 10 and 45000 in steps of 10, where OrdMag(*x*) is the number of barcodes
#' with UMI count above a threshold 10-fold less than the 99th percentile of the top
#' *x* barcodes.
#'
#' @param counts Integer vector of total UMI counts per barcode.
#'
#' @returns Integer scalar *x*.
#'
#' @export
estimate.expected.cells <- function(counts) {
  df <- data.frame(
    x = seq_along(counts),
    logcounts = sort(log10(counts), decreasing = TRUE),
    stat = NA
  ) %>%
    filter(x <= 45000)

  for (i in seq.int(from = 10, to = nrow(df), by = 10)) {
    threshold <- .ordmag(df$logcounts[seq_len(i)])
    df$stat[i] <- (sum(df$logcounts > threshold) - df$x[i])^2 / df$x[i]
  }

  return(df$x[match(min(df$stat, na.rm = TRUE), df$stat)])
}


#' k-means cell calling algorithm
#'
#' @description
#' Modified implementation of CellRanger ARC cell calling algorithm.
#'
#' Algorithm steps:
#' 1. Compute log10-transformed total counts for each modality per barcode
#' 2. Filter to barcodes with count >= 1 for every modality
#' 3. Deduplicate barcodes with identical counts across all modalities
#' 4. Determine initial thresholds for cell-containing barcodes using OrdMag algorithm
#'    for each modality and compute centroids for cell-containing and empty barcodes
#' 5. Peform k-means clustering (k = 2) initialised using centroids from OrdMag-based
#'    thresholds and update barcode labels based on cluster assignment
#' 6. Project barcode labels from deduplicated set to full set
#'
#' Key differences from CellRanger ARC implementation are:
#' * Lack of initial ATAC prefiltering steps
#' * Generalisation to more than 2 modalities
#'
#' @param matrix.list List of raw count matrices.
#' @param n.expected.cells Integer scalar. Targeted cell recovery.
#' @param ordmag.quantile Numeric scalar (default 0.99). Quantile used for OrdMag function.
#' @param ordmag.ratio Numeric scalar (default 10). Ratio used for OrdMag function.
#' @param seed Integer scalar (default 100). Random seed.
#' @param verbose Logical scalar (default `TRUE`). Print progress messages.
#'
#' @returns Character vector of cell-containing barcodes.
#'
#' @importFrom dplyr filter everything select distinct mutate if_else if_all group_by summarise across starts_with arrange case_match left_join pull
#' @importFrom tibble rownames_to_column column_to_rownames
#'
#' @export
kmeans.cell.caller <- function(matrix.list,
                               n.expected.cells,
                               ordmag.quantile = 0.99,
                               ordmag.ratio = 10,
                               seed = 100,
                               verbose = TRUE) {
  if (!is.list(matrix.list)) stop("'matrix.list' must be a list")
  if (is.null(names(matrix.list))) names(matrix.list) <- as.character(seq_along(matrix.list))

  if (verbose) message("Running k-means cell calling algorithm using ", paste(toupper(names(matrix.list)), collapse = ", "))

  set.seed(as.integer(seed))

  matrix.list <- lapply(matrix.list, function(x) x[, colSums(x) > 0])

  logcounts <- Map(
    f = function(x, i) {
      data.frame(log10(colSums(x))) %>%
        setNames(paste("logcounts", i, sep = ".")) %>%
        tibble::rownames_to_column("barcode")
    },
    x = matrix.list,
    i = names(matrix.list)
  ) %>%
    Reduce(
      x = .,
      f = function(x, y) dplyr::full_join(x, y, by = "barcode"),
      init = data.frame(barcode = rownames(matrix.list[[1]]))
    ) %>%
    dplyr::filter(dplyr::if_all(dplyr::everything(), function(x) !is.na(x)))

  deduplicated <- logcounts %>%
    dplyr::select(dplyr::starts_with("logcounts")) %>%
    dplyr::distinct()
  centers <- deduplicated %>%
    dplyr::mutate(
      cell = dplyr::if_else(
        dplyr::if_all(
          dplyr::everything(),
          function(x) {
            threshold <- .ordmag(
              x = sort(x, decreasing = TRUE)[seq_len(n.expected.cells)],
              q = ordmag.quantile,
              r = ordmag.ratio
            )
            return(x > threshold)
          }
        ),
        TRUE,
        FALSE
      )
    ) %>%
    dplyr::group_by(cell) %>%
    dplyr::summarise(dplyr::across(dplyr::starts_with("logcounts"), mean)) %>%
    dplyr::arrange(cell) %>%
    tibble::column_to_rownames("cell") %>%
    as.matrix()
  deduplicated <- deduplicated %>%
    dplyr::mutate(
      cell = dplyr::case_match(kmeans(x = ., centers = centers)$cluster, 1 ~ FALSE, 2 ~ TRUE)
    )

  cells <- dplyr::left_join(
    x = logcounts,
    y = deduplicated,
    by = paste("logcounts", names(matrix.list), sep = ".")
  ) %>%
    dplyr::filter(cell) %>%
    dplyr::pull("barcode")

  return(cells)
}


#' Wrapper function to run EmptyDrops algorithm on a single modality
#'
#' @param x Raw count matrix.
#' @param modality Character scalar. Modality name.
#' @param n.expected.cells Integer scalar. Targeted cell recovery.
#' @param lower Integer scalar. Upperbound of total counts for barcodes to use for estimating ambient profile.
#' @param ignore Integer scalar. Upperbound of total counts for barcodes to consider as non-empty droplets.
#' @param niters Integer scalar (default 10000). Number of Monte Carlo iterations.
#' @param seed Integer scalar (default 100). Random seed.
#' @param verbose Logical scalar (default `TRUE`). Print progress messages.
#' @param ... Additional arguments passed to `DropletUtils::testEmptyDrops`.
#'
#' @returns Data frame with columns:
#' * PValue: p-value for each barcode.
#' * FDR: false discovery rate for each barcode.
#' * Limited: logical vector indicating whether p-values are limited by the number of Monte Carlo
#' iterations.
#'
#' @importFrom DropletUtils testEmptyDrops
#'
#' @export
.emptydrops <- function(x,
                        modality,
                        lower,
                        ignore,
                        niters = 10000,
                        seed = 100,
                        verbose = TRUE,
                        ...) {
  if (verbose) message("Running EmptyDrops on ", toupper(modality), " with ", niters, " iterations")
  set.seed(as.integer(seed))
  res <- DropletUtils::testEmptyDrops(
    x,
    lower = lower,
    ignore = ignore,
    niters = niters,
    ...
  )
  res$FDR <- p.adjust(res$PValue, method = "BH")
  rownames(res) <- colnames(x)
  return(res)
}


#' Multimodal EmptyDrops cell calling algorithm
#'
#' @description
#' Modified implementation of CellRanger variant of EmptyDrops cell calling algorithm for multimodal
#' data.
#'
#' Algorithm steps:
#' 1. Prefilter barcodes with high total counts for ALL modalities using OrdMag-derived thresholds
#' based on the number of expected cells (these are considered to definitely be cell-containing)
#' 2. Determine 'lower' and 'ignore' thresholds for each modality based on user-defined parameters
#' `ind.min`, `umi.min`, `umi.min.frac.median` and `cand.max.n`
#' 3. Run EmptyDrops on each modality separately, excluding prefiltered barcodes and barcodes with ranks
#' higher than the upperbound rank used for estimating the ambient profiles `ind.max` (the latter are
#' considered to only contain counts due to index hopped reads and therefore are not used for estimating
#' the ambient profile)
#' 4. Check if p-values for any modality are limited by the number of Monte Carlo iterations
#' and rerun EmptyDrops with an increased number of iterations if necessary
#' 5. Add default values for p-values and FDRs for barcodes excluded from EmptyDrops
#' 6. Combine p-values across modalities using Fisher method
#' 7. Adjust combined p-values using Benjamini-Hochberg method
#'
#' @param matrix.list List of raw count matrices.
#' @param n.expected.cells Integer scalar. Targeted cell recovery.
#' @param ordmag.quantile Numeric scalar (default 0.99). Quantile used for OrdMag function.
#' @param ordmag.ratio Numeric scalar (default 10). Ratio used for OrdMag function.
#' @param umi.min Integer scalar (default 500). Minimum UMI count per barcode.
#' @param umi.min.frac.median Numeric scalar (default 0.01). Minimum fraction of median UMI count
#' per barcode.
#' @param ind.min Integer scalar (default 45000). Lowerbound rank of barcodes to use for estimating
#' ambient profile.
#' @param ind.max Integer scalar (default 90000). Upperbound rank of barcodes to use for estimating
#' ambient profile.
#' @param cand.max.n Integer scalar (default 20000). Maximum number of barcodes to consider as
#' potentially cell-containing.
#' @param niters Integer scalar (default 10000). Number of Monte Carlo iterations.
#' @param max.attempts Integer scalar (default 3). Maximum number of attempts to rerun EmptyDrops if
#' p-values are limited.
#' @param seed Integer scalar (default 100). Random seed.
#' @param verbose Logical scalar (default `TRUE`). Print progress messages.
#' @param ... Additional arguments passed to `DropletUtils::testEmptyDrops`.
#'
#' @returns Data frame with columns:
#' * PValue.modality: p-value for each modality.
#' * FDR.modality: false discovery rate for each modality.
#' * PValue: combined p-value.
#' * FDR: combined false discovery rate.
#' Additionally, the data frame has an attribute 'limited' which is a logical vector indicating
#' whether p-values for each modality are limited.
#'
#' @importFrom dplyr if_else if_any starts_with mutate select
#'
#' @export
emptydrops.multimodal <- function(matrix.list,
                                  n.expected.cells,
                                  ordmag.quantile = 0.99,
                                  ordmag.ratio = 10,
                                  umi.min = 500,
                                  umi.min.frac.median = 0.01,
                                  cand.max.n = 20000,
                                  ind.min = 45000,
                                  ind.max = 90000,
                                  niters = 10000,
                                  max.attempts = 3,
                                  seed = 100,
                                  verbose = TRUE,
                                  ...) {
  if (!is.list(matrix.list)) stop("'matrix.list' must be a list")
  if (is.null(names(matrix.list))) names(matrix.list) <- as.character(seq_along(matrix.list))

  # store original barcode order
  barcodes <- colnames(matrix.list[[1]])

  # prefilter barcodes with high counts in all modalities which are high likely to be cells
  # regardless of whether their feature profile differs significantly from the ambient profile
  totals <- lapply(matrix.list, Matrix::colSums)
  prefiltered <- lapply(
    totals,
    function(total, n.expected.cells) {
      threshold <- .ordmag(
        x = sort(log10(total), decreasing = TRUE)[seq_len(n.expected.cells)],
        q = ordmag.quantile,
        r = ordmag.ratio
      )
      return(names(total)[log10(total) >= threshold])
    },
    n.expected.cells = n.expected.cells
  ) %>%
    Reduce(intersect, ., init = barcodes)
  # rank barcodes by total counts and determine for each modality:
  # - lower: the count threshold below which barcodes are used to estimate the ambient profile
  # - ignore: the count threshold below which barcodes are not considered as potentially cell-containing
  # (but not necessarily used to estimate the ambient profile)
  # - exclude: the count threshold below which barcodes are assumed to only contain index hopped reads
  ranked <- lapply(totals, sort, decreasing = TRUE)
  lower <- lapply(ranked, `[`, ind.min)
  ignore <- lapply(ranked, function(ranked) {
    max(umi.min, round(umi.min.frac.median * median(ranked[prefiltered])), ranked[cand.max.n])
  })
  exclude <- lapply(ranked, function(ranked) {
    names(ranked[seq.int(from = ind.max, to = length(ranked))])
  })
  # remove barcodes with total counts below ignore threshold from prefiltered set
  prefiltered <- mapply(
    function(total, minimum) {
      names(total)[total >= minimum]
    },
    total = totals,
    minimum = ignore,
    SIMPLIFY = FALSE
  ) %>%
    Reduce(intersect, ., init = prefiltered)
  # remove prefiltered and excluded barcodes from count matrices prior to running EmptyDrops
  filtered <- mapply(function(x, exclude) {
    x[, !colnames(x) %in% c(exclude, prefiltered)]
  }, x = matrix.list, exclude = exclude, SIMPLIFY = FALSE)

  # initialise results lists
  ed.res <- vector("list", length(filtered)) %>% setNames(names(filtered))
  is.limited <- rep(TRUE, length(filtered)) %>% setNames(names(filtered))
  attempt <- 1
  # run EmptyDrops on each modality which increasing number of iterations until all modalities are
  # no longer limited or maximum number of attempts is reached
  while (any(is.limited) && attempt <= max.attempts) {
    if (attempt > 1 && verbose) message("The following modalities require more iterations: ", paste(toupper(names(ed.res)[is.limited]), collapse = ", "))
    ed.res[is.limited] <- mapply(
      .emptydrops,
      x = filtered[is.limited],
      modality = names(filtered[is.limited]),
      lower = lower[is.limited],
      ignore = ignore[is.limited],
      MoreArgs = list(niters = attempt * niters, seed = seed, verbose = verbose, ...),
      SIMPLIFY = FALSE
    )
    is.limited <- sapply(ed.res, function(x) any(x$Limited[which(x$FDR >= 0.001)]))
    attempt <- attempt + 1
  }
  if (any(is.limited)) message("Some modalities are still limited; consider increasing 'niters' and/or 'max.attempts'.")

  # add default values for prefiltered and excluded barcodes to results
  # all barcodes in this set are assigned PValue = NA (not used in EmptyDrops)
  # prefiltered barcodes are assigned FDR = 0 (definitely cell-containing)
  # excluded barcodes are assigned FDR = NA (definitely not cell-containing)
  default.res <- lapply(
    exclude,
    function(exclude, prefiltered) {
      data.frame(
        PValue = rep(NA, length(c(exclude, prefiltered))),
        FDR = c(rep(NA, length(exclude)), rep(0, length(prefiltered))),
        row.names = c(exclude, prefiltered)
      )
    },
    prefiltered = prefiltered
  )
  res <- mapply(function(x, y) {
    df <- rbind(x[c("PValue", "FDR")], y[c("PValue", "FDR")])
    return(df[barcodes, ]) # return results in original barcode order
  }, x = ed.res, y = default.res, SIMPLIFY = FALSE)

  if (verbose && length(res) > 1) message("Aggregating results across modalities")
  # combine PValue and FDR from each modality and aggregate results
  pvalue <- lapply(res, `[[`, "PValue") %>%
    setNames(paste("PValue", names(.), sep = ".")) %>%
    as.data.frame()
  fdr <- lapply(res, `[[`, "FDR") %>%
    setNames(paste("FDR", names(.), sep = ".")) %>%
    as.data.frame()
  out <- cbind(pvalue, fdr) %>%
    dplyr::mutate(
      PValue = apply(dplyr::select(., dplyr::starts_with("PValue.")), MARGIN = 1, FUN = function(x) ifelse(any(is.na(x)), NA, ifelse(length(x) > 1, metap::sumlog(x)$p, x))),
      FDR = p.adjust(PValue, method = "BH")
    ) %>%
    dplyr::mutate(
      FDR = dplyr::if_else(dplyr::if_any(dplyr::starts_with("FDR."), function(x) x == 0), 0, FDR)
    )
  # store limited flag for each modality as a data frame attribute
  attr(out, "limited") <- is.limited
  return(out)
}


#' Combine results of multiple QC filters
#'
#' Takes a list or data frame of QC filtering results (logical vectors with
#' 1 element per cell barcode) and combines the results to identify barcodes which
#' have passed all filters.
#'
#' @param filters List or data frame of logical vectors.
#' @param missing Logical scalar (default FALSE). Value to convert NA values
#' in `filters` to prior to combining results.
#'
#' @returns Logical vector.
#'
#' @importFrom magrittr or
#'
#' @export
combine.filters <- function(filters, missing = FALSE) {
  if ((!is.list(filters) || !is.data.frame(filters)) && !all(sapply(filters, is.logical))) stop("filters must be a list of logical vectors")

  for (i in seq_along(filters)) {
    filters[[i]] <- ifelse(is.na(filters[[i]]), missing, filters[[i]])
  }

  return(Reduce(magrittr::or, filters, init = FALSE))
}


#' Generate scatter plot
#'
#' Makes a scatter plot.
#'
#' @param data Data frame or object coercible to data frame.
#' @param x Column in `data` to determine x position.
#' @param y Column in `data` to determine y position.
#' @param colour Column in `data` to determine colour.
#' @param ... Fixed aesthetics to pass to [ggplot2::geom_point()].
#' @param title Character scalar specifying plot title.
#' @param x.lab Character scalar specifying x-axis label.
#' @param y.lab Character scalar specifying y-axis label.
#' @param colour.lab Character scalar specifying colour legend label.
#' @param log Log-transform axis scales:
#'  * `TRUE` to transform all axes
#'  * Character vector with values 'x' and/or 'y' to transform specified axes
#'
#' @returns Plot object.
#'
#' @importFrom ggplot2 ggplot aes geom_point labs theme_minimal scale_x_log10 scale_y_log10
#' @importFrom scales label_number cut_short_scale
#'
#' @export
scatter.plot <- function(data,
                         x,
                         y,
                         colour,
                         ...,
                         title = NULL,
                         x.lab = NULL,
                         y.lab = NULL,
                         colour.lab = NULL,
                         log = NULL) {
  data <- data %>%
    as.data.frame()

  plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = {{ x }}, y = {{ y }}, colour = {{ colour }})) +
    ggplot2::geom_point(...) +
    ggplot2::labs(title = title, x = x.lab, y = y.lab, colour = colour.lab) +
    ggplot2::theme_minimal()
  if (!is.null(log)) {
    if (is.logical(log) && log) log <- c("x", "y")
    if ("x" %in% log) plot <- plot + ggplot2::scale_x_log10(labels = scales::label_number(scale_cut = append(scales::cut_short_scale(), 1, 1)))
    if ("y" %in% log) plot <- plot + ggplot2::scale_y_log10(labels = scales::label_number(scale_cut = append(scales::cut_short_scale(), 1, 1)))
  }

  return(plot)
}


#' Generate violin/beeswarm plot
#'
#' Makes a combined violin plot and beeswarm plot.
#'
#' @inheritParams scatter.plot
#' @param ... Fixed aesthetics to pass to [ggplot2::geom_quasirandom()].
#' @param points Logical scalar (default `TRUE`). Add beeswarm points.
#'
#' @returns Plot object.
#'
#' @importFrom ggplot2 ggplot aes geom_violin labs theme_minimal scale_y_log10
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom scales label_number cut_short_scale
#'
#' @export
violin.plot <- function(data,
                        x,
                        y,
                        colour,
                        ...,
                        points = TRUE,
                        title = NULL,
                        x.lab = NULL,
                        y.lab = NULL,
                        colour.lab = NULL,
                        log = NULL) {
  data <- data %>%
    as.data.frame()

  plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = {{ x }}, y = {{ y }}, colour = {{ colour }})) +
    ggplot2::geom_violin(colour = "grey50") +
    ggplot2::labs(title = title, x = x.lab, y = y.lab, colour = colour.lab) +
    ggplot2::theme_minimal()
  if (points) {
    plot <- plot +
      ggbeeswarm::geom_quasirandom(...)
  }
  if (!is.null(log)) {
    if (is.logical(log) && log) log <- "y"
    if ("x" %in% log) stop("Unable to log transform categorical x-axis scale")
    if ("y" %in% log) plot <- plot + ggplot2::scale_y_log10(labels = scales::label_number(scale_cut = append(scales::cut_short_scale(), 1, 1)))
  }

  return(plot)
}


#' Generate histogram plot
#'
#' Makes a histogram plot.
#'
#' @param data Data frame or object coercible to data frame.
#' @param x Column in `data` to determine x position.
#' @param y Function to determine y position after [ggplot2::stat_bin()]; see [ggplot2::after_stat()].
#' @param fill Column in `data` to determine fill.
#' @param ... Fixed aesthetics to pass to [ggplot2::geom_histogram()].
#' @param title Character scalar specifying plot title.
#' @param x.lab Character scalar specifying x-axis label.
#' @param y.lab Character scalar specifying y-axis label.
#' @param fill.lab Character scalar specifying fill legend label.
#' @param log Log-transform axis scales:
#' * `TRUE` to transform all axes
#' * Character vector with values 'x' and/or 'y' to transform specified axes
#'
#' @returns Plot object.
#'
#' @importFrom ggplot2 ggplot geom_histogram labs theme_minimal scale_x_log10 scale_y_log10 after_stat
#' @importFrom scales label_number cut_short_scale
#'
#' @export
hist.plot <- function(data,
                      x,
                      y = ggplot2::after_stat(count),
                      fill,
                      ...,
                      title = NULL,
                      x.lab = NULL,
                      y.lab = NULL,
                      fill.lab = NULL,
                      log = FALSE) {
  plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = {{ x }}, y = {{ y }}, fill = {{ fill }})) +
    ggplot2::geom_histogram(position = "identity", ...) +
    ggplot2::labs(title = title, x = x.lab, y = y.lab, fill = fill.lab) +
    ggplot2::theme_minimal()
  if (!is.null(log)) {
    if (is.logical(log) && log) log <- c("x", "y")
    if ("x" %in% log) plot <- plot + ggplot2::scale_x_log10(labels = scales::label_number(scale_cut = append(scales::cut_short_scale(), 1, 1)))
    if ("y" %in% log) plot <- plot + ggplot2::scale_y_log10(labels = scales::label_number(scale_cut = append(scales::cut_short_scale(), 1, 1)))
  }

  return(plot)
}


#' Generate bar plot
#'
#' Makes a bar plot.
#'
#' @param data Data frame or object coercible to data frame.
#' @param pos Column in `data` to determine bar position.
#' @param fill Column in `data` to determine bar fill.
#' @param ... Fixed aesthetics to pass to [ggplot2::geom_col()].
#' @param title Character scalar specifying plot title.
#' @param pos.lab Character scalar specifying categorical axis label.
#' @param fill.lab Character scalar specifying fill legend label.
#' @param log Log-transform axis scales:
#'  * `TRUE` to transform all axes
#'  * Character vector with values 'x' and/or 'y' to transform specified axes
#'
#' @returns Plot object.
#'
#' @importFrom dplyr group_by summarise n
#' @importFrom ggplot2 ggplot geom_col scale_y_discrete labs theme_minimal scale_x_log10
#' @importFrom scales label_number cut_short_scale
#'
#' @export
bar.plot <- function(data,
                     pos,
                     fill = pos,
                     ...,
                     title = NULL,
                     pos.lab = NULL,
                     fill.lab = NULL,
                     log = NULL) {
  data <- data %>%
    as.data.frame() %>%
    dplyr::group_by({{ pos }}, {{ fill }}) %>%
    dplyr::summarise(count = dplyr::n())

  plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = count, y = {{ pos }}, fill = {{ fill }})) +
    ggplot2::geom_col(...) +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::labs(title = title, y = pos.lab, fill = fill.lab) +
    ggplot2::theme_minimal()
  if (!is.null(log)) {
    if (is.logical(log) && log) log <- "x"
    if ("x" %in% log) plot <- plot + ggplot2::scale_x_log10(labels = scales::label_number(scale_cut = append(scales::cut_short_scale(), 1, 1)))
    if ("y" %in% log) stop("Unable to log transform categorical y-axis scale")
  }

  return(plot)
}


#' Generate barcode rank plot
#'
#' Makes a barcode rank plot coloured by fraction of barcodes at each rank called
#' as cells using the output of [DropletUtils::barcodeRanks()].
#'
#' @param bcrank.res DataFrame output of [DropletUtils::barcodeRanks()].
#' @param cells Character vector of cell-containing barcodes.
#' @param ... Arguments to pass to [scatter.plot()].
#'
#' @returns Plot object.
#'
#' @export
bcrank.plot <- function(bcrank.res, cells, ...) {
  data <- bcrank.res %>%
    as.data.frame() %>%
    mutate(cell = if_else(rownames(.) %in% cells, TRUE, FALSE)) %>%
    group_by(rank, total) %>%
    summarise(fraction.cells = sum(cell) / n()) %>%
    arrange(rank, desc(fraction.cells))

  plot <- scatter.plot(
    data,
    x = rank,
    y = total,
    colour = fraction.cells,
    log = TRUE,
    ...
  )

  return(plot)
}


#' Generate feature log-counts histogram plots
#'
#' Makes log-counts histogram plots for selected features.
#'
#' @param matrix Matrix or matrix-like object containing counts with
#'  features as row names.
#' @param features.pattern Character vector of regular expressions to select
#'  features for plotting.
#' #' @param ... Arguments to pass to [hist.plot()].
#'
#' @returns Plot object.
#'
#' @importFrom dplyr bind_cols everything
#' @importFrom tidyr pivot_longer
#'
#' @export
logcounts.histplot <- function(matrix,
                               features.pattern = "^[a-zA-Z]+_",
                               ...) {
  data <- matrix[grepl(pattern = features.pattern, x = rownames(matrix)), ] %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    tidyr::pivot_longer(cols = dplyr::everything(), names_to = "Feature", values_to = "count")

  plot <- hist.plot(
    data = data,
    x = count,
    fill = Feature,
    log = "x",
    ...
  )

  return(plot)
}


#' Generate feature log-counts ridge plots
#'
#' Makes log-counts distribution plots for selected features, optionally
#' grouping by sample.
#'
#' @param matrix Matrix or matrix-like object containing counts with
#'  features as row names.
#' @param ... Fixed aesthetics to pass to [ggridges::geom_density_ridges()].
#' @param features.pattern Character vector of regular expressions to select
#'  features for plotting.
#' @param pseudocount Numeric scalar (default 1). Added to counts before log
#'  transformation (to avoid infinite values).
#' @param x.lower Numeric scalar (default 0.3). Lower limit of x-axis scale
#'  (choose a value greater than 0 to avoid loss of resolution due to inflated
#'  counts close to zero; set to `NA` to use ggplot2 default).
#' @param x.upper Numeric scalar (default maximum log-count in `matrix`).
#'  Upper limit of x-axis scale (set to `NA` for to use ggplot2 default).
#' @param samples Character vector of the same length as number of
#'  cells specifying sample of origin. If specified, will group by sample.
#' @param title Character scalar specifying plot title.
#' @param x.lab Character scalar specifying x-axis label.
#' @param y.lab Character scalar specifying y-axis label.
#'
#' @returns Plot object.
#'
#' @importFrom dplyr bind_cols mutate across everything
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 aes geom_density_ridges scale_y_discrete xlim labs theme_minimal theme
#' @importFrom ggridges geom_density_ridges
#'
#' @export
logcounts.ridgeplot <- function(matrix,
                                ...,
                                features.pattern = "^[a-zA-Z]+_",
                                pseudocount = 1,
                                x.lower = 0.3,
                                x.upper = max(log10(matrix + pseudocount)),
                                samples = NULL,
                                title = NULL,
                                x.lab = sprintf("log10(UMI count + %d)", pseudocount),
                                y.lab = ifelse(is.null(samples), "Feature", "Sample")) {
  data <- matrix %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ log10(.x + pseudocount)))
  if (!is.null(samples)) data <- dplyr::bind_cols(data, Sample = samples)
  data <- tidyr::pivot_longer(data, cols = dplyr::matches(features.pattern) & !dplyr::matches("^Sample$"), names_to = "Feature", values_to = "count")

  plot <- ggplot2::ggplot(data = data)
  if (is.null(samples)) {
    plot <- plot + ggridges::geom_density_ridges(mapping = ggplot2::aes(x = count, y = Feature), ...)
  } else {
    plot <- plot +
      ggridges::geom_density_ridges(mapping = ggplot2::aes(x = count, y = Sample, fill = Sample), ...) +
      ggplot2::facet_wrap(~Feature, scales = "free_x", ncol = 2, labeller = "label_both")
  }
  plot <- plot +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::xlim(x.lower, x.upper) +
    ggplot2::labs(
      title = title,
      x = x.lab,
      y = y.lab
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(strip.text = ggplot2::element_text(hjust = 0))

  return(plot)
}


#' Generate demultiplexing scatter biplots
#'
#' Makes scatter biplots of log-counts for selected features, coloured by demultiplexing classification.
#'
#' @param matrix Matrix or matrix-like object containing counts with
#' features as row names.
#' @param res A data frame with one row per cell barcode and character columns 'HTO' and 'Classification'
#' indicating HTO assignment(s) (comma-separated if multiple) and final classification respectively.
#' @param feature.1 Character scalar specifying first feature to plot.
#' @param feature.2 Character scalar specifying second feature to plot.
#' @param ... Fixed aesthetics to pass to [scatter.plot()].
#'
#' @returns Plot object.
#'
#' @importFrom dplyr arrange bind_cols
#' @importFrom stringr str_split
#' @importFrom rlang sym
#'
#' @export
demux.plot <- function(matrix,
                       res,
                       feature.1,
                       feature.2,
                       ...) {
  data <- matrix %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    dplyr::bind_cols(res)
  keep <- sapply(
    stringr::str_split(res$HTO, ","),
    function(x) {
      if (length(x) == 1) {
        out <- any(c(feature.1, feature.2, "negative", "uncertain") %in% x)
      } else {
        out <- all(c(feature.1, feature.2) %in% x)
      }
      return(out)
    }
  )
  data <- data[keep, ] %>%
    dplyr::arrange(Classification)

  plot <- scatter.plot(
    data = data,
    x = !!rlang::sym(feature.1),
    y = !!rlang::sym(feature.2),
    colour = Classification,
    x.lab = paste(feature.1, "UMI count"),
    y.lab = paste(feature.2, "UMI count"),
    colour.lab = "Classification",
    log = TRUE,
    ...
  )

  return(plot)
}


#' Generate UMAP plot of known and predicted multiplets
#'
#' Makes UMAP plot of known and predicted multiplets. If multiple dimensionality reductions
#' are present in the input Seurat object, a weighted nearest neighbour UMAP projection is
#' computed.
#'
#' @param x A Seurat object contianing precomputed dimensionality reduction (PCA/LSI).
#' @param res A data frame with one row per cell barcode and boolean columns 'known' and 'predicted'
#' indicating whether the barcode is a known or predicted multiplet respectively.
#' @param ... Fixed aesthetics to pass to [ggplot2::geom_point()].
#'
#' @returns Plot object.
#'
#' @importFrom SeuratObject AddMetaData Reductions FetchData
#' @importFrom Seurat FindMultiModalNeighbors FindNeighbors RunUMAP
#' @importFrom dplyr case_when
#' @importFrom ggplot2 ggplot aes geom_point scale_colour_manual labs theme_minimal
#'
#' @export
multiplet.plot <- function(x,
                           res,
                           ...) {
  if (!is(x, "Seurat")) stop("x must be an object of class 'Seurat'")

  x <- SeuratObject::AddMetaData(x, res)
  x$Type <- dplyr::case_when(
    x$known ~ "known",
    x$predicted ~ "predicted",
    TRUE ~ "singlet"
  )

  if (length(SeuratObject::Reductions(x)) > 1) {
    x <- Seurat::FindMultiModalNeighbors(
      x,
      reduction.list = SeuratObject::Reductions(x),
      dims.list = lapply(SeuratObject::Reductions(x), function(x) seq.int(from = if (x == "lsi") 2 else 1, to = 30)),
      weighted.nn.name = "nn",
      verbose = FALSE
    )
  } else {
    x <- Seurat::FindNeighbors(
      x,
      reduction = SeuratObject::Reductions(x),
      dims = seq.int(from = if (SeuratObject::Reductions(x) == "lsi") 2 else 1, to = 30),
      return.neighbor = TRUE,
      graph.name = "nn",
      verbose = FALSE
    )
  }
  x <- Seurat::RunUMAP(x, nn.name = "nn", spread = 10, verbose = FALSE) # increase spread to improve multiplet separation for visualisation

  plot <- ggplot2::ggplot(data = SeuratObject::FetchData(x, vars = c("umap_1", "umap_2", "Type")), mapping = ggplot2::aes(x = umap_1, y = umap_2, colour = Type)) +
    ggplot2::geom_point(...) +
    ggplot2::scale_colour_manual(values = c("singlet" = "grey80", "known" = "#d82526", "predicted" = "#ffc156")) +
    ggplot2::labs(title = "Multiplet assignment", x = "UMAP1", y = "UMAP2") +
    ggplot2::theme_minimal()

  return(plot)
}


#' Generate TSS enrichment plot
#'
#' Makes a plot of mean enrichment scores by distance from TSS for a set of cell barcodes.
#'
#' @param x Object of class 'Seurat' or 'ChromatinAssay' containing TSS enrichment scores.
#' Alternatively, a sparse matrix of TSS enrichment scores with cell barcodes as rows and distance
#' from TSS as columns.
#' @param assay Character scalar specifying assay to use in 'Seurat' object. (default is to
#' use the default assay).
#' @param cells Character vector of cell barcodes (default is to use all barcodes).
#' @param ... Fixed aesthetics to pass to [ggplot2::geom_line()].
#'
#' @returns Plot object.
#'
#' @importFrom Seurat DefaultAssay GetAssayData
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot aes geom_line labs theme_minimal
#'
#' @export
tss.plot <- function(x,
                     assay = NULL,
                     cells = NULL,
                     ...) {
  if (is(x, "Seurat")) {
    assay <- assay %||% Seurat::DefaultAssay(x)
    x <- x[[assay]]
    if (!is(x, "ChromatinAssay")) stop("Assay must be an object of class 'ChromatinAssay'")
  }
  if (is(x, "ChromatinAssay")) {
    x <- Seurat::GetAssayData(x, layer = "positionEnrichment")[["TSS"]]
    if (is.null(x)) stop("No TSS enrichment scores found in assay")
  }
  if (!is(x, "sparseMatrix")) {
    stop("x must be an object of class 'Seurat' or 'ChromatinAssay' containing TSS enrichment scores or a sparse matrix of TSS enrichment scores")
  }
  if (!is.null(cells)) x <- x[cells, ]

  data <- Matrix::colMeans(x)
  data <- data %>%
    data.frame(mean.enrichment = .) %>%
    tibble::rownames_to_column("distance") %>%
    dplyr::mutate(distance = as.integer(distance))
  plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = distance, y = mean.enrichment)) +
    ggplot2::geom_line(...) +
    ggplot2::labs(title = "Tn5 insertion frequency", x = "Distance from TSS (bp)", y = "Relative frequency") +
    ggplot2::theme_minimal()

  return(plot)
}


#' Generate histogram of fragment size distribution
#'
#' Makes a histogram of fragment size distribution for a set of cell barcodes.
#' Uses a representative region of the genome to avoid reading all fragments into
#' memory.
#'
#' @param fragments Fragments object.
#' @param cells Character vector of cell barcodes (default is to use all barcodes in fragments object).
#' @param region Character vector specifying genomic regions to use (default is to use first million bases
#' from chromosomes 1-22 and X).
#' @param ... Fixed aesthetics to pass to [ggplot2::geom_freqpoly()].
#'
#' @returns Plot object.
#'
#' @importFrom Signac GetFragmentData GetReadsInRegion
#' @importFrom ggplot2 ggplot aes geom_freqpoly scale_y_continuous xlim labs theme_minimal
#'
#' @export
fragments.plot <- function(fragments,
                           cells = NULL,
                           region = paste(sprintf("chr%s", c(as.character(1:22), "X")), "1", "1000000", sep = "-"),
                           ...) {
  if (!is(fragments, "Fragment")) stop("fragments must be an object of class 'Fragment'")

  data <- Signac:::GetReadsInRegion(
    cellmap = Signac::GetFragmentData(fragments, slot = "cells"),
    region = region,
    tabix.file = Signac::GetFragmentData(fragments, slot = "path"),
    cells = cells
  )
  plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = length, y = ggplot2::after_stat(density), text = paste("percent: ", after_stat(density) * 100))) +
    ggplot2::geom_freqpoly(...) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::xlim(20, 700) +
    ggplot2::labs(title = "Fragment length distribution", x = "Length (bp)", y = "Fragments (%)") +
    ggplot2::theme_minimal()

  return(plot)
}
