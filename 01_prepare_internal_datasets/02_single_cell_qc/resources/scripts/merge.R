#!/bin/env -S Rscript --vanilla


# ==============================
# COMMAND LINE OPTIONS
# ==============================
# Define options
DOC <- "
Merge multimodal count matrices into a single 10x Genomics HDF5 file

Usage:
  merge.R [--gex=<path>] [--atac=<path>] [--adt=<path>] [--adt-prefix=<prefix>] --output=<path> [options]

Arguments:
  REQUIRED
  -o --output=<path>          Path to output matrix (10x Genomics HDF5 format)

  OPTIONAL:
  --gex=<path>                Path to GEX matrix (10x Genomics Matrix Market format)
  --atac=<path>               Path to ATAC matrix (10x Genomics Matrix Market format)
  --adt=<path>                Path to ADT matrix (10x Genomics Matrix Market format or BarCounter CSV format)
  --adt-prefix=<prefix>       Prefix for ADT features [default: ADT_]

Options:
  -h --help                   Show this screen
"

# Parse options
opt <- docopt::docopt(DOC)

# Logging options
logger::log_layout(logger::layout_glue)
logger::log_warnings()
logger::log_errors()


# ==============================
# SETUP
# ==============================
logger::log_info("Initialising")

suppressPackageStartupMessages({
  library(Matrix)
})

source("utils.R")

# Check input
if (is.null(opt$gex) && is.null(opt$atac) && is.null(opt$adt)) {
  stop("At least one of --gex, --atac or --adt must be provided")
}


# ==============================
# SCRIPT
# ==============================
counts <- list()
if (!is.null(opt$gex)) {
  logger::log_info("Loading GEX matrix: {opt$gex}")
  counts[["Gene Expression"]] <- get.10x.matrix(file = opt$gex, type = "Gene Expression")
  counts[["Gene Expression"]] <- counts[["Gene Expression"]][, colSums(counts[["Gene Expression"]]) > 0] # remove cell barcodes with zero counts
}
if (!is.null(opt$atac)) {
  logger::log_info("Loading ATAC matrix: {opt$atac}")
  counts[["Peaks"]] <- get.10x.matrix(file = opt$atac, type = "Peaks")
  counts[["Peaks"]] <- counts[["Peaks"]][, colSums(counts[["Peaks"]]) > 0] # remove cell barcodes with zero counts
}
if (!is.null(opt$adt)) {
  logger::log_info("Loading ADT matrix: {opt$adt}")
  if (grepl(pattern = ".mtx.gz$", x = opt$adt)) {
    counts[["Antibody Capture"]] <- get.10x.matrix(file = opt$adt, type = "Antibody Capture")
    counts[["Antibody Capture"]] <- counts[["Antibody Capture"]][grep(pattern = paste0("^", opt$adt_prefix), rownames(counts[["Antibody Capture"]])), ]
  } else {
    counts[["Antibody Capture"]] <- get.barcounter.matrix(file = opt$adt, features.pattern = paste0("^", opt$adt_prefix))
  }
  counts[["Antibody Capture"]] <- counts[["Antibody Capture"]][, colSums(counts[["Antibody Capture"]]) > 0] # remove cell barcodes with zero counts
}
counts <- counts[sapply(counts, function(x) nrow(x) > 0)] # remove empty matrices

logger::log_info("Merging count matrices")
bc.common <- Reduce(intersect, lapply(counts, colnames))
counts <- lapply(counts, function(x) x[, bc.common])

logger::log_info("Saving output")
if (!dir.exists(dirname(opt$output))) dir.create(dirname(opt$output), recursive = TRUE)
feature.counts <- sapply(counts, nrow)
DropletUtils::write10xCounts(
  path = opt$output,
  x = Reduce(rbind, counts),
  gene.type = rep(names(feature.counts), times = feature.counts),
  genome = "GRCh38",
  version = "3",
  overwrite = TRUE
) |>
  suppressWarnings() |>
  suppressMessages()
logger::log_success("Output path: {opt$output}")
