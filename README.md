# DPhil Thesis Reproducibility Repository

This repository contains the analysis code and workflow for my DPhil thesis:

**Characterising the effect of therapeutic intervention to discover pathogenic immune interactions in multiple sclerosis**

## Overview

This repository provides the code, configuration files and software environment specifications necessary to reproduce the analyses presented in the thesis. The workflow is organized into four main stages:
1. Raw data preprocessing, quality control, normalisation and integration of CITE-seq datasets generated for this study
2. Preprocessing, normalisation and integration of published external datasets used for validation analyses
3. Clustering and cell type annotation
4. Downstream analyses

## Repository Structure

```
ocrelizumab_paper/
├── 01_prepare_internal_datasets/   # Snakemake pipelines for preprocessing, QC, normalisation and integration of internal datasets
├── 02_prepare_external_datasets/   # Notebooks for preprocessing, normalisation and integration of external datasets
├── 03_annotate/                    # Notebooks for clustering and cell type annotation
├── 04_analyse/                     # Notebooks for downstream analyses
├── data/                           # Data directory
├── envs/                           # Conda environment specifications
├── results/                        # Analysis outputs
├── scripts/                        # Utility scripts
└── slurm/                          # SLURM job submission scripts
```

## Prerequisites

### System Requirements
- Linux operating system (tested on Ubuntu 22.04)
- High-performance computing cluster with NVIDIA CUDA-compatible GPU recommended (tested on [MRC WIMM JADE HPC cluster](https://www.imm.ox.ac.uk/facilities/ccb-high-performance-computing/jade-hpc-cluster))

### Software Dependencies

#### Snakemake Pipelines
The Snakemake pipelines in [`01_prepare_internal_datasets/`](01_prepare_internal_datasets/) contain their own Conda environment specifications and installation instructions.

#### Global Environment
All software dependencies for the Jupyter notebooks are specified in [`envs/environment.yaml`](envs/environment.yaml)

To create the Conda environment, first install Anaconda or Miniconda (see [installation instructions](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html)), then run:
```bash
conda env create -f envs/environment.yaml
conda activate ocrelizumab_paper
```

Next, to install the Jupyter kernel for running R notebooks, run:
```bash
Rscript -e 'IRkernel::installspec(name = "ir-ocrelizumab_paper", displayname = "R 4.3 (ocrelizumab_paper)")'
```

## Data Availability

### Input Data
Internal datasets:
- FASTQ files for internal CITE-seq datasets will be deposited in EGA (accession pending)
- Processed data for internal CITE-seq datasets have been deposited in GEO (accession GSE316688):
    - Unfiltered per-sample GEX and ADT count matrices
    - QC-filtered and annotated datasets (Seurat/RDS and MuData/H5MU formats) - download to `data/processed/cite_seq/cohort_treatment_naive/annotated/` and `data/processed/cite_seq/cohort_nonresponders/annotated/`

External datasets:
- Cantoni et al. (2025): Synapse (accession syn51730532) - download to `data/processed/external/cantoni/source/`
- Absinta et al. (2021): GEO (accession GSE180759) - download to `data/processed/external/lesion_rims/absinta/source/`
- Lerma-Martin et al. (2024): EGA (accession EGAC50000000231) - download to `data/processed/external/lesion_rims/lerma_martin/source/`
- Kaufmann et al. (2021): GEO (accession GSE144744) - download to `data/processed/external/kaufmann/source/`

### Metadata

Internal datasets:
- Anonymised donor-level metadata for internal CITE-seq datasets are provided in [`data/metadata/`](data/metadata/).

External datasets:
- Additional manually-curated metadata for external datasets (from supplementary materials of original publications) are provided in [`data/processed/external/`](data/processed/external/).

### Reference Data

Specific reference data files used in the analyses are provided in [`data/ref/`](data/ref/). Additional reference files (e.g. public datasets and CellTypist models for reference mapping and cell type annotation as specified in notebooks) should be downloaded to this directory.

### cNMF Gene Programs and starCAT Spectra

Precomputed cNMF gene programs and starCAT spectra from the discovery cohort dataset are provided in [`results/tables/cite_seq/cohort_treatment_naive/gep/`](results/tables/cite_seq/cohort_treatment_naive/gep/).

## General Notes

- All Snakemake pipelines and notebooks are numbered and should be executed in order; subsequent stages depend on outputs from previous stages.
- Computationally intensive analysis steps were run on an HPC cluster using SLURM job submission scripts provided in the [`slurm/`](slurm/) directory; these steps are indicated within notebooks at the relevant points. These shell scripts will need to be adapted to the specific HPC environment and job scheduler in use.
- Relative file paths are used in all notebooks; ensure that the repository structure is maintained when running the code. Shell scripts include a `BASEDIR` variable that should be set to the root directory of the copy of repository on the local system.

## License

All code in this repository is published under the MIT License (see the [LICENSE](LICENSE) file for details).

## Contact

For questions about the code, please contact Redwan Farooq ([redwan.farooq@ndcn.ox.ac.uk](mailto:redwan.farooq@ndcn.ox.ac.uk)).

---

**Last updated**: 24/04/2026