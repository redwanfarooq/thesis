# Stage 1: Raw data preprocessing, quality control, normalisation and integration of CITE-seq datasets generated for this study

This stage consists of sequential Snakemake pipelines for processing internal CITE-seq datasets. The pipelines are self-contained with their own Conda environment specifications and installation instructions.

The datasets were generated with three experimental protocols:
- `cohort_01`: Chromium Next GEM Single Cell 3' v3.1 with TotalSeq-A antibodies
- `cohort_02`: Chromium Next GEM Single Cell 3' v3.1 with TotalSeq-B antibodies
- `cohort_03`: Chromium GEM-X Single Cell 3' v4 with TotalSeq-B antibodies

`cohort_01` and `cohort_02` comprised the discovery cohort (`cohort_treatment_naive`) and `cohort_03` comprised the non-responder validation cohort (`cohort_nonresponders`).

Accordingly, configuration files (`config`) and metadata files (`metadata`) are provided separately for each cohort within each pipeline directory. Copy the appropriate configuration and metadata files to the main `config` and `metadata` directories within each pipeline before running. Paths to input and output files/directories should be modified as necessary.

The pipelines are configured to run on a SLURM-based HPC cluster; however, this may need to be adapted to the specific HPC environment and job scheduler in use (`profile/config.yaml`).

1. [`01_preprocessing/`](01_preprocessing/): 
   - Demultiplexing BCL files with bcl2fastq
   - Alignment and quantification with STARsolo (GEX) and BarCounter (ADT/HTO)
   - Per-sample sequencing and alignment QC reports FastQC and MultiQC

2. [`02_single_cell_qc/`](02_single_cell_qc/):
   - Ambient decontamination with CellBender
   - Cell calling with EmptyDrops
   - Cell hashing demultiplexing with demuxmix
   - Doublet detection with scDblFinder
   - Automated outlier-based QC filtering with scater
   - Per-sample cell-level QC reports with custom scripts
   
3. [`03_single_cell_multi/`](03_single_cell_multi/):
   - Merge per-sample datasets into a combined Seurat object
   - Gene symbol standardisation to HGNC-approved symbols
   - Normalisation (RNA: LogNormalize, ADT: CLR)
   - HVG selection
   - Regression of technical covariates (percentage mitochondrial genes, cell cycle scores)
   - PCA
   - Integration with Harmony

**Notes:**

1. [`01_preprocessing/`](01_preprocessing/) was originally run using BCL files as input. However, raw data files for this project have been deposited in EGA as FASTQ files. To run the pipeline from FASTQ files, provide a modified runs summary table in accordance with the [pipeline documentation](01_preprocessing/README.md). Note the specific requirements for FASTQ file naming and directory structure.
2. [`02_single_cell_qc/`](02_single_cell_qc/) can be run using the unfiltered per-sample count matrices provided in GEO (accession GSE316688) as input, bypassing the need to run [`01_preprocessing/`](01_preprocessing/). Note the files downloaded from GEO will need to be reorganised into the expected file naming and directory structure (i.e. one directory per sample containing the GEX and ADT count matrices with matrices renamed to `matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz` and `Tag_Counts.csv` with no prefixes).