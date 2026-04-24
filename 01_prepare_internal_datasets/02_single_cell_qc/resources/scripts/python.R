# ==============================
# Python wrapper functions
# ==============================


require(reticulate)


# Python wrapper function for running COMPOSITE algorithm
reticulate::source_python("run_composite.py")


#' Run COMPOSITE
#'
#' This function is an R wrapper to run the COMPOSITE algorithm for joint multiplet calling on count matrices
#' from one or more modalities.
#'
#' @param ... One or more sparse matrices containing count data from different modalities. The matrices must be
#' named according to the modality with the following options: `"gex"`, `"atac"`, or `"adt"`.
#' @param fragments A character string specifying the path to a CellRanger-formatted single cell ATAC fragments file.
#' Required if an ATAC count matrix is provided. Default is `NULL`.
#' @param tmpdir A character string specifying the directory to write temporary files. Default is `tempdir()`.
#' @param verbose A logical value indicating whether to print progress messages. Default is `TRUE`.
#'
#' @return A data frame with the following columns:
#' - `classification`: A character vector indicating the classification of each cell barcode as either "singlet" or "multiplet".
#' - `probability`: A numeric vector indicating the posterior probability of each cell barcode being a multiplet.
#'
#' @importFrom reticulate import
#' @importFrom Matrix writeMM
#' @importFrom dplyr case_match
#'
#' @export
run.composite <- function(...,
                          fragments = NULL,
                          tmpdir = tempdir(),
                          verbose = TRUE) {
  x <- list(...)
  if (length(x) == 0) stop("Invalid input: no count matrices provided")
  if (!all(sapply(x, function(x) is(x, "sparseMatrix")))) stop("Invalid input: all count matrices must be sparse matrices")
  if (length(x) > 1 && !Reduce(identical, lapply(x, colnames))) stop("Invalid input: all count matrices must have the same cell barcodes")
  names(x) <- sapply(names(x), function(x) match.arg(tolower(x), choices = c("gex", "atac", "adt")))
  if (any(is.null(names(x)))) stop("Invalid count matrix type(s): '...' arguments must be named 'gex', 'atac', or 'adt'")
  if ("atac" %in% names(x) && is.null(fragments)) stop("Invalid input: 'fragments' must be provided with ATAC count matrix")
  if (!is.null(fragments) && !file.exists(fragments)) stop("Invalid input: 'fragments' file does not exist")
  names(x) <- dplyr::case_match(
    names(x),
    "gex" ~ "RNA",
    "atac" ~ "ATAC",
    "adt" ~ "ADT"
  )

  # store cell barcodes
  barcodes <- colnames(x[[1]])

  # COMPOSITE does not work well with ATAC peak count matrix so convert to gene activity matrix
  if ("ATAC" %in% names(x)) {
    if (verbose) message("Computing gene activity matrix for ATAC")
    annotation <- Signac::GetGRangesFromEnsDb(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86) |> suppressWarnings()
    GenomeInfoDb::seqlevelsStyle(annotation) <- "UCSC"
    GenomeInfoDb::genome(annotation) <- "hg38"
    obj <- SeuratObject::CreateSeuratObject(
      counts = Signac::CreateChromatinAssay(
        counts = x[["ATAC"]],
        sep = c(":", "-"),
        fragments = fragments,
        annotation = annotation,
        verbose = FALSE
      ),
      assay = "ATAC"
    )
    x[["ATAC"]] <- Signac::GeneActivity(obj, assay = "ATAC", gene.id = TRUE, verbose = FALSE)
  }

  if (verbose) message("Writing count matrices to temporary files")
  x <- lapply(
    x,
    function(x) {
      file <- tempfile(tmpdir = tmpdir, fileext = ".mtx")
      Matrix::writeMM(x, file)
      return(file)
    }
  )

  if (verbose) message("Running COMPOSITE with the following modalities: ", paste(names(x), collapse = ", "))
  x$file <- tempfile(tmpdir = tmpdir, fileext = ".tsv")
  # run COMPOSITE and parse stdout to extract goodness-of-fit scores for each modality
  res <- do.call(run_composite, x) |>
    reticulate::py_capture_output() |>
    stringr::str_split_1(pattern = "\n") |>
    stringr::str_match(pattern = "The (?<modality>[A-Z]+) modality goodness-of-fit score is: (?<score>[0-9.]+)") |>
    na.omit()
  scores <- as.numeric(x = res[, "score"]) %>% setNames(res[, "modality"])

  out <- read.delim(x$file, sep = "\t", header = FALSE, row.names = barcodes, col.names = c("classification", "probability")) |>
    dplyr::mutate(classification = ifelse(classification == 1, "multiplet", "singlet"))
  attr(out, "scores") <- scores

  return(out)
}
