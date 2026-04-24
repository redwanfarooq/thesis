# Description
Snakemake pipeline for preprocessing single cell sequencing data
- Modularised workflow can be modified and/or extended for different library preparation protocols
- Add as a submodule in a bioinformatics project GitHub repository
```
git submodule add https://github.com/redwanfarooq/preprocessing preprocessing
```
- Update submodule to the latest version
```
git submodule update --remote preprocessing
```

# Required software
1. Global environment
    - [Snakemake >=v7.31](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
    - [docopt >=v0.6](https://github.com/docopt/docopt)
    - [pandas >=v2.0](https://pandas.pydata.org/docs/getting_started/install.html)
    - [loguru >=v0.7](https://github.com/Delgan/loguru)
    - [h5py >=3.12](https://docs.h5py.org/en/latest/build.html)
2. Specific modules
    - [bcl2fastq >=v2.20](https://sapac.support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html)
    - [Seqtk >= 1.3](https://github.com/lh3/seqtk)
    - [Cell Ranger >=v9.0](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation)
    - [Cell Ranger ATAC >=v2.0](https://software.10xgenomics.com/single-cell-atac/software/pipelines/latest/installation)
    - [Cell Ranger ARC >=v2.0](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/installation)
    - [STAR >=v2.7.11a](https://github.com/alexdobin/STAR)
    - [samtools >=v1.17](http://www.htslib.org)
    - [chromap >=v0.2.5](https://github.com/haowenz/chromap)
    - [htslib >=v1.18](http://www.htslib.org)
    - [MACS2 >=v2.2.9](https://github.com/macs3-project/MACS/wiki/Install-macs2)
    - [BarCounter](https://github.com/AllenInstitute/BarCounter-release)
    - [FastQC >=v0.11](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
    - [MultiQC >=v1.14](https://multiqc.info/docs/getting_started/installation/)
    - [quarto >=v1.4](https://quarto.org/docs/get-started/)
    - [panel >=v1.3.8](https://panel.holoviz.org/getting_started/installation.html)
    - [R >=v4.3](https://cran.r-project.org)
        * [Signac v1.10.0](https://CRAN.R-project.org/package=Signac)

# Setup
1. Install software for global environment (requires Anaconda or Miniconda - see [installation instructions](https://conda.io/projects/conda/en/stable/user-guide/install/index.html))
    - Download [environment YAML](/resources/envs/snakemake.yaml)
    - Create new conda environment from YAML
    ```
    conda env create -f snakemake.yaml
    ```
2. Install software for specific module(s)
    - Manually install required software from source and check that executables are available in **PATH** (using `which`) *and/or*
    - Create new conda environments with required software from YAML (as above - download [environment YAMLs](/resources/envs)) *and/or*
    - Check that required software is available to load as environment modules (using `module avail`)
3. Set up pipeline configuration file **config/config.yaml** (see comments in file for detailed instructions)
4. Set up profile configuration file **profile/config.yaml** (see comments in file for detailed instructions)

# Run
1. Activate global environment
```
conda activate snakemake
```
2. Execute **run.py** in root directory

# Input
Pipeline requires the following input files/folders:

## General

**REQUIRED:**

1. Folder(s) containing input BCL files
- One folder for each Illumina sequencing run with standard directory structure (**RunInfo.xml** must be present at top level) 
- Folders should ideally be named according to default convention for the system e.g. **YYYYMMDD_InstrumentID_RunNumber_FlowCellID**, but any folder naming convention ending in an underscore followed by a unique ID will suffice

*or*

Folder(s) containing input FASTQ files
- One folder for each Illumina sequencing run with subfolders containing FASTQ files from each library type 
- Folders should ideally be named according to default convention for the system e.g. **YYYYMMDD_InstrumentID_RunNumber_FlowCellID**, but any folder naming ending in an underscore followed by a unique ID will suffice
- Subfolders should be named according to library type (must match __exactly__ with **lib_type** field entry in runs summary table)
- FASTQ files should be named according to default convention e.g. **SampleID_Sx_Lxxx_Rx_001.fastq.gz**
2. Runs summary table in delimited file format (e.g. TSV, CSV) with the following required fields (with headers):
- **run**: run folder name
- **format**: file type (options: BCL, FASTQ)
- **lib_type**: library type (options: GEX, ATAC, ADT, HTO, CRISPR, BCR, TCR)
- **sample_id**: sample ID
- **sample_index**: *either* index name *or* literal i7 index sequence - only required if **format** is BCL
- **lane**: *either* lane number(s) (separated with spaces if more than one lane used) *or* * (for all lanes) - only required if **format** is BCL
- **sample_index2**: literal i5 index sequence (if applicable) - only required if **format** is BCL, dual indexing used and **sample_index** is literal i7 index sequence

**OPTIONAL:**

3. Index kit CSV files with the following required fields (with headers) - must be provided if any index names are used in **sample_index** field of runs summary table:
- **index_name**: index set name (e.g. SI-TT-A1); must be first field
- **\***: 1 or more fields specifying index sequences in the following order:
    - *Dual index kits*: i7 sequence, i5 sequence (forward strand), i5 sequence (reverse complement)
    - *Single index kits*: i7 sequence (additional fields if more than 1 sequence per index set)

## Module-specific

### gex: GEX only protocol

**REQUIRED:**

1. Reference files:
- STAR genome reference package
- Cell barcode whitelist

### atac: ATAC only protocol

**REQUIRED:**

1. Reference files:
- chromap genome reference and index
- Cell barcode whitelist

### gex_atac: 10X multiome (GEX + ATAC) protocol

**REQUIRED:**

1. Reference files:
- STAR genome reference package
- chromap genome reference and index
- Cell barcode whitelist (GEX)
- Cell barcode whitelist (ATAC)

### cite_seq: CITE-seq protocol

**REQUIRED:**

1. Reference files:
- STAR genome reference package
- Cell barcode whitelist (GEX)
2. Antibody tag list in CSV format with the following required fields (without headers):
- Tag sequence (length 15nt) - must begin at first base in read 2 (if leading bases are present, FASTQ files must be trimmed e.g. TotalSeq-B and TotalSeq-C antibodies) 
- Tag name

### tea_seq: TEA-seq/DOGMA-seq protocol

**REQUIRED:**

1. Reference files:
- STAR genome reference package
- chromap genome reference and index
- Cell barcode whitelist (GEX)
- Cell barcode whitelist (ATAC)
2. Antibody tag list in CSV format with the following required fields (without headers):
- Tag sequence (length 15nt) - must begin at first base in read 2 (if leading bases are present, FASTQ files must be trimmed e.g. TotalSeq-B and TotalSeq-C antibodies) 
- Tag name

### gex_fb_cellranger: 10X GEX +/- feature barcoding protocol

**REQUIRED:**

1. Reference files:
- Cell Ranger genome reference package
2. Feature reference in CSV format - see [specifications](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis#feature-ref) (if using feature barcoding)

### atac_cellranger: 10X ATAC only protocol

**REQUIRED:**

1. Reference files:
- Cell Ranger genome reference package

### gex_atac_cellranger: 10X multiome (GEX + ATAC) protocol

**REQUIRED:**

1. Reference files:
- Cell Ranger genome reference package

### vdj_cellranger: 10X 5' immune profiling protocol

**REQUIRED:**
1. Reference files:
- Cell Ranger VDJ reference package

### tea_seq_cellranger: TEA-seq/DOGMA-seq protocol

**REQUIRED:**

1. Reference files:
- Cell Ranger genome reference package
- Cell barcode whitelist (GEX)
2. Antibody tag list in CSV format with the following required fields (without headers):
- Tag sequence (length 15nt) - must begin at first base in read 2 (if leading bases are present, FASTQ files must be trimmed e.g. TotalSeq-B and TotalSeq-C antibodies)
- Tag name

### gex_fb_vdj_cellranger: 10X GEX +/- feature barcoding +/- immune profiling protocol (with optional sample hashing)

**REQUIRED:**
1. Reference files:
- Cell Ranger genome reference package
- Cell Ranger VDJ reference package (if using immune profiling)
2. Feature reference in CSV format - see [specifications](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis#feature-ref) (if using feature barcoding)
3. Sample hashing table in in delimited file format (e.g. TSV, CSV) with the following required fields (with headers) (if using sample hashing):
- **sample_id**: sample ID (must match **sample_id** field in runs summary table)
- **hash_id**: hashed sample ID (must be unique for each hashed sample)
- **ocm_barcode_ids/hashtag_ids/cmo_ids**: OCM barcode/hashtag/CMO IDs (if antibody hashtags are used, must match **id** field in feature reference CSV) 

# Output
Output directory will be created in specified location with subfolders containing the output of each software tool specified in the module.

# Modules

## Available modules
- gex
- atac
- gex_atac
- cite_seq
- tea_seq
- gex_fb_cellranger
- atac_cellranger
- gex_atac_cellranger
- vdj_cellranger
- tea_seq_cellranger
- gex_fb_vdj_cellranger

## Adding new module
1. Add entry to module rule specifications file **config/modules.yaml** with module name and list of rule names
2. Add additional rule definition files in **modules/rules** folder (if needed)
- Rule definition file **must** also assign a list of pipeline target files generated by the rule to a variable with the same name as the rule
- Rule definition file **must** have the same file name as the rule with the file extension **.smk**
3. Execute **run.py** in root directory with `--update` flag (needs to be repeated if there are any further changes to the module rule specification in **config/modules.yaml**)