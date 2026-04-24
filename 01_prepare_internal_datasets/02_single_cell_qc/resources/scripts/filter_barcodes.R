#!/bin/env -S Rscript --vanilla


# ==============================
# COMMAND LINE OPTIONS
# ==============================
# Define options
DOC <- "
Filter merged multimodal count matrices by cell barcode 

Usage:
  filter_barcodes.R --input=<path> --barcodes=<path> --output=<path> [options]

Arguments:
  REQUIRED
  -i --input=<path>           Path to input matrix (10x Genomics HDF5 format)
  -b --barcodes=<path>        Path to text file containing cell barcodes to keep (one per line)
  -o --output=<path>          Path to output matrix (10x Genomics HDF5 format)

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

# ==============================
# SCRIPT
# ==============================
logger::log_info("Loading barcodes: {opt[['--barcodes']]}")
valid <- vroom::vroom(
  opt[["--barcodes"]],
  delim = "\t",
  col_names = FALSE,
  show_col_types = FALSE
) |>
  dplyr::pull(1)

logger::log_info("Loading count matrices: {opt[['--input']]}")
counts <- get.10x.h5(opt[["--input"]])
logger::log_info("Filtering count matrices")
counts <- lapply(counts, function(x, valid) x[, valid], valid = valid)

logger::log_info("Saving output")
if (!dir.exists(dirname(opt[["--output"]]))) dir.create(dirname(opt[["--output"]]), recursive = TRUE)
feature.counts <- sapply(counts, nrow)
DropletUtils::write10xCounts(
  path = opt[["--output"]],
  x = Reduce(rbind, counts),
  gene.type = rep(names(feature.counts), times = feature.counts),
  genome = "GRCh38",
  version = "3",
  overwrite = TRUE
) |>
  suppressWarnings() |>
  suppressMessages()
logger::log_success("Output path: {opt[['--output']]}")
