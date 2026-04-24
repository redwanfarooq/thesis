# Stage 2: Preprocessing, normalisation and integration of published external datasets used for validation analyses

This stage consists of sequential notebooks for processing published external validation datasets. The notebooks rely on the global Conda environment specification provided in [`envs/environment.yaml`](../envs/environment.yaml).

The same general preprocessing steps were applied to all datasets where applicable:
- Gene symbol standardisation to HGNC-approved symbols
- Normalisation (RNA: LogNormalize, ADT: CLR)
- HVG selection
- Regression of technical covariates (percentage mitochondrial genes, cell cycle scores)
- PCA
- Integration with Harmony

1. [`01_gene_translation.ipynb`](01_gene_translation.ipynb) - gene symbol standardisation
2. [`02_preprocess_cantoni.ipynb`](02_preprocess_cantoni.ipynb) - preprocess Cantoni et al. (2025) dataset
3. [`03_preprocess_lesion_rims.ipynb`](03_preprocess_lesion_rims.ipynb) - preprocess Absinta et al. (2021) and Lerma-Martin et al. (2024) lesion rim datasets
4. [`04_preprocess_kaufmann.ipynb`](04_preprocess_kaufmann.ipynb) - preprocess Kaufmann et al. (2021) dataset
5. [`05_convert_seurat.ipynb`](05_convert_seurat.ipynb) - convert processed AnnData/MuData objects to Seurat objects for interoperability