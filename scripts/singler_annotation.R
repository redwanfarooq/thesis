#!/bin/env -S Rscript --vanilla


# ==============================
# COMMAND LINE OPTIONS
# ==============================
# Define options
DOC <- "
SingleR annotation

Usage:
  singler_annotation.R --query=<query> --reference=<reference> --vars=<var[;var...]> --output=<output> [options]

Arguments:
  REQUIRED
  --query=<query>                         Path to RDS or QS file containing Seurat object
  --reference=<reference>                 Reference dataset in celldex package
  --vars=<var;[var...]>                   Semicolon-separated list of metadata variables to transfer from reference
  --output=<output>                       Path to the output TSV file

  OPTIONAL
  -t --threads=<int>                      Number of threads [default: 1]
  -l --log=<file>                         Path to log file [default: singler_annotation.log]
  --query-assay=<assay>                   Query assay name [default: RNA]
  --query-layer=<layer>                   Query layer name [default: data]

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
  library(Seurat)
  library(BPCells)
  library(celldex)
  library(SingleR)
})

if (!dir.exists(dirname(params$output))) dir.create(dirname(params$output), recursive = TRUE)


# ==============================
# SCRIPT
# ==============================
logger::log_info("Loading query dataset")
# Load query dataset
if (grepl(pattern = "\\.qs$", x = params$query, ignore.case = TRUE)) {
  seu <- qs::qread(file = params$query, nthreads = as.integer(params$threads))
} else {
  seu <- readRDS(file = params$query)
}
sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(logcounts = LayerData(seu, assay = params$query_assay, layer = params$query_layer))
)

logger::log_info("Loading reference dataset")
# Load reference dataset
if (!params$reference %in% ls("package:celldex")) {
  stop("Invalid reference dataset: ", params$reference, ". Available datasets in celldex package: ", paste(ls("package:celldex"), collapse = ", "))
}
reference <- get(params$reference, envir = asNamespace("celldex"))()

# Run SingleR annotation
logger::log_info("Transferring annotations from reference dataset: {paste(params$vars, collapse = ', ')}")
annotation <- set_names(params$vars) %>%
  map(
    .f = function(x, test = sce, ref = reference, nthreads = as.integer(params$threads), verbose = !params$quiet) {
      if (verbose) message("Predicting ", x)
      res <- SingleR::SingleR(
        test = test,
        ref = ref,
        labels = ref[[x]],
        restrict = intersect(rownames(test), rownames(ref)),
        BPPARAM = BiocParallel::MulticoreParam(workers = nthreads)
      )
      labels <- res[, "pruned.labels", drop = FALSE] %>%
        setNames(x) %>%
        as.data.frame()
      return(labels)
    }
  ) %>%
  do.call(cbind, .)

# Save output
logger::log_info("Saving output")
write.table(
  annotation,
  file = params$output,
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = TRUE
)
logger::log_success("Output path: {params$output}")
