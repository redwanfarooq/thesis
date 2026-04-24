#!/bin/env -S Rscript --vanilla


# ==============================
# COMMAND LINE OPTIONS
# ==============================
# Define options
DOC <- "
Reference mapping for multisample single cell experiments

Usage:
  reference_map.R --query=<query> --reference=<reference> --vars=<var[;var...]> --output=<output> [options]

Arguments:
  REQUIRED
  --query=<query>                         Path to RDS or QS file containing merged Seurat object with multimodal count matrices and metadata
  --reference=<reference>                 Path to RDS or QS file containing reference Seurat object with multimodal count matrices and metadata
  --vars=<var;[var...]>                   Semicolon-separated list of metadata variables to transfer from reference
  --output=<output>                       Path to the output TSV file

  OPTIONAL
  -t --threads=<int>                      Number of threads [default: 1]
  -l --log=<file>                         Path to log file [default: reference_map.log]
  --exp=<name>                            Name of metadata field to use for experiment variable [default: orig.ident]
  --query-assay=<assay>                   Query assay name [default: RNA]
  --reference-assay=<assay>               Reference assay name [default: RNA]
  --reference-reduction=<reduction>       Reference dimensionality reduction name [default: spca]
  --normalisation-method=<method>         Normalisation method used for reference assay ('LogNormalize' or 'SCT') [default: LogNormalize]
  --k=<int>                               Number of nearest neighbours to use for reference mapping [default: 50]
  --n-dims=<int>                          Number of dimensions to use for reference mapping [default: 50]

Options:
  -h --help                               Show this screen
  -q --quiet                              Do not print logging messages to console
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
params <- lapply(opt, function(x) if (is.character(x)) stringr::str_split_1(string = x, pattern = ";") else x)

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
  library(dplyr)
  library(purrr)
  library(future)
  library(Seurat)
  library(BPCells)
})

options(future.globals.maxSize = 50000 * 1024^2)

if (!dir.exists(dirname(params$output))) dir.create(dirname(params$output), recursive = TRUE)


# ==============================
# SCRIPT
# ==============================
plan(multicore, workers = as.integer(params$threads))

logger::log_info("Loading query dataset")
params$normalisation_method <- params$normalisation_method %>% match.arg(choices = c("LogNormalize", "SCT"))
# Load query dataset
if (grepl(pattern = "\\.qs$", x = params$query, ignore.case = TRUE)) {
  seu <- qs::qread(file = params$query, nthreads = as.integer(params$threads))
} else {
  seu <- readRDS(file = params$query)
}
logger::log_info("Loading reference dataset")
# Load reference dataset
if (grepl(pattern = "\\.qs$", x = params$reference, ignore.case = TRUE)) {
  ref <- qs::qread(file = params$reference, nthreads = as.integer(params$threads))
} else {
  ref <- readRDS(file = params$reference)
}
if (params$normalisation_method == "SCT" && Seurat:::IsSCT(seu[[params$query_assay]])) {
  umi.assay <- slot(object = seu[[params$query_assay]], name = "SCTModel.list") %>%
    map(.f = function(x) slot(object = x, name = "umi.assay")) %>%
    unique() %>%
    unlist()
  params$query_assay <- c(umi.assay, params$query_assay)
}
seu <- DietSeurat(object = seu, assays = params$query_assay)
ref <- DietSeurat(object = ref, assays = params$reference_assay, dimreducs = params$reference_reduction)

# Get experiment variable and split query dataset by experiment
params$exp <- params$exp %>% match.arg(choices = colnames(seu[[]]))
seu <- SplitObject(object = seu, split.by = params$exp)
# Get metadata variables to transfer
params$vars <- params$vars %>%
  match.arg(choices = colnames(ref[[]]), several.ok = TRUE) %>%
  set_names()


# Precompute nearest neighbours graph
logger::log_info("Computing nearest neighbours graph")
ref <- FindNeighbors(
  object = ref,
  reduction = params$reference_reduction,
  dims = seq_len(as.integer(params$n_dims)),
  graph.name = "annoy.neighbors",
  k.param = as.integer(params$k),
  cache.index = TRUE,
  return.neighbor = TRUE,
  l2.norm = TRUE
)
# Compute anchors
logger::log_info("Computing anchors between reference and query datasets")
anchors <- map(
  .x = seu,
  .f = function(query,
                reference = ref,
                normalization.method = params$normalisation_method,
                reference.reduction = params$reference_reduction,
                dims = seq_len(as.integer(params$n_dims)),
                verbose = !params$quiet) {
    FindTransferAnchors(
      query = query,
      reference = reference,
      normalization.method = normalization.method,
      reference.reduction = reference.reduction,
      reference.neighbors = "annoy.neighbors",
      dims = dims,
      verbose = verbose
    )
  }
)
# Transfer data
logger::log_info("Transferring data from reference dataset: {paste(params$vars, collapse = ', ')}")
data <- map(
  .x = anchors,
  .f = function(anchorset,
                reference = ref,
                vars = params$vars,
                verbose = !params$quiet) {
    TransferData(
      anchorset = anchorset,
      reference = reference,
      refdata = as.list(vars),
      verbose = verbose
    ) %>%
      map2(
        .y = params$vars,
        .f = function(x, var) x[, "predicted.id", drop = FALSE] %>% set_names(var)
      ) %>%
      bind_cols()
  }
) %>%
  bind_rows()


# Save output
logger::log_info("Saving output")
write.table(
  data,
  file = params$output,
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = TRUE
)
logger::log_success("Output path: {params$output}")
