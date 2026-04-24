#!/bin/env -S Rscript --vanilla


# ==============================
# COMMAND LINE OPTIONS
# ==============================
# Define options
DOC <- "
Generate peak count matrix

Usage:
  count_peaks.R --fragments=<fragments> --peaks=<peaks> [--threads=<threads>] [options]

Arguments:
  -f --fragments=<fragments>  Path to fragments TSV file (required)
  -p --peaks=<peaks>          Path to peaks BED file (required)
  -t --threads=<threads>      Number of threads [default: 1]

Options:
  -h --help                   Show this screen
"

# Parse options
opt <- docopt::docopt(DOC)
outdir <- file.path(dirname(opt[["--fragments"]]), "raw_feature_bc_matrix")

# Logging options
logger::log_layout(logger::layout_glue)
logger::log_warnings()
logger::log_errors()


# ==============================
# SETUP
# ==============================
logger::log_info("Initialising")

suppressPackageStartupMessages({
  library(future)
  library(dplyr)
  library(Signac)
  library(GenomicRanges)
  library(Matrix)
})

logger::log_info("Running with {as.integer(opt[['--threads']])} threads")
plan(multicore, workers = as.integer(opt[["--threads"]]))


# ==============================
# SCRIPT
# ==============================
logger::log_info("Loading fragments: {opt[['--fragments']]}")
fragments <- CreateFragmentObject(opt[["--fragments"]])

logger::log_info("Loading peaks: {opt[['--peaks']]}")
features <- read.table(
  file = opt[["--peaks"]],
  header = FALSE,
  sep = "\t"
) %>%
  select(1:3) %>%
  transmute(
    chr = V1,
    start = V2,
    stop = V3,
    id = sprintf("%s:%d-%d", chr, start, stop),
    symbol = sprintf("%s:%d-%d", chr, start, stop),
    type = "Peaks"
  ) %>%
  select(id, symbol, type, chr, start, stop) %>%
  distinct(id, .keep_all = TRUE)
peaks <- makeGRangesFromDataFrame(features, starts.in.df.are.0based = TRUE)

logger::log_info("Counting fragments in peaks per cell")
peaks.mat <- FeatureMatrix(fragments, peaks)

logger::log_info("Saving output")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
write.table(
  features,
  file = file.path(outdir, "features.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  data.frame(barcodes = colnames(peaks.mat)),
  file = file.path(outdir, "barcodes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

writeMM(
  peaks.mat,
  file = file.path(outdir, "matrix.mtx")
) %>% invisible()

logger::log_success("Output path: {outdir}")
