# Stage 3: Clustering and cell type annotation

This stage consists of sequential notebooks for clustering and cell type annotation of all datasets. The notebooks rely on the global Conda environment specification provided in [`envs/environment.yaml`](../envs/environment.yaml).

1. [`01_cluster_annotate_cohort_treatment_naive.ipynb`](01_cluster_annotate_cohort_treatment_naive.ipynb) - semi-supervised consensus clustering and annotation of the discovery cohort dataset using scTriangulate
2. [`02_subcluster_annotate_cohort_treatment_naive.ipynb`](02_subcluster_annotate_cohort_treatment_naive.ipynb) - fine-grained subclustering and annotation of the discovery cohort dataset
3. [`03_build_reference.ipynb`](03_build_reference.ipynb) - construction of a PBMC reference from the discovery cohort dataset for annotation of other datasets
4. [`04_subcluster_annotate_cohort_nonresponders.ipynb`](04_subcluster_annotate_cohort_nonresponders.ipynb) - fine-grained subclustering and annotation of the non-responder cohort dataset based on reference mapping
5. [`05_subcluster_annotate_cantoni.ipynb`](05_subcluster_annotate_cantoni.ipynb) - fine-grained subclustering and annotation of the CSF/blood validation dataset based on reference mapping
6. [`06_subcluster_annotate_lesion_rims.ipynb`](06_subcluster_annotate_lesion_rims.ipynb) - fine-grained subclustering and annotation of the brain tissue validation dataset based on reference mapping
7. [`07_subcluster_annotate_kaufmann.ipynb`](07_subcluster_annotate_kaufmann.ipynb) - fine-grained subclustering and annotation of the natalizumab validation dataset based on reference mapping
8. [`08_cluster_names.ipynb`](08_cluster_names.ipynb) - standardised cell type nomenclature
9. [`09_save_annotated.ipynb`](09_save_annotated.ipynb) - export annotated and filtered Seurat and AnnData objects for use in downstream analyses