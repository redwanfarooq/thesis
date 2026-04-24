# Stage 4: Downstream Analyses

This stage consists of sequential notebooks for downstream analyses of annotated datasets. The notebooks rely on the global Conda environment specification provided in [`envs/environment.yaml`](../envs/environment.yaml).

1. [`01_meld.ipynb`](01_meld.ipynb) - MELD analysis to compute perturbation scores
2. [`02_cnmf.ipynb`](02_cnmf.ipynb) - cNMF analysis to infer cell type-specific gene programs and compute program scores
3. [`03_gep_visualisation.ipynb`](03_gep_visualisation.ipynb) - visualisation of program score distributions
4. [`04_gep_pathway_enrichment.ipynb`](04_gep_pathway_enrichment.ipynb) - pathway enrichment analysis of top 100 genes per program
5. [`05_rf_shap.ipynb`](05_rf_shap.ipynb) - random forest regression of perturbation scores vs program scores and compute SHAP feature importance
6. [`06_shap_program_perturbation.ipynb`](06_shap_program_perturbation.ipynb) - visualisation of SHAP values vs program scores vs perturbation scores
7. [`07_totalvi.ipynb`](07_totalvi.ipynb) - totalVI denoising and differential expression analysis
8. [`08_clinical_demographics.ipynb`](08_clinical_demographics.ipynb) - visualisation and summary of clinical demographics
9. [`09_cohort_treatment_naive_global.ipynb`](09_cohort_treatment_naive_global.ipynb) - discovery cohort global analysis
10. [`10_cohort_treatment_naive_b_cells.ipynb`](10_cohort_treatment_naive_b_cells.ipynb) - discovery cohort B cell analysis
11. [`11_cohort_treatment_naive_t_cells.ipynb`](11_cohort_treatment_naive_t_cells.ipynb) - discovery cohort T cell analysis
12. [`12_cohort_treatment_naive_cd20dim_t_cells.ipynb`](12_cohort_treatment_naive_cd20dim_t_cells.ipynb) - discovery cohort CD20dim T cells analysis
13. [`11_cohort_treatment_naive_t_cells.ipynb`](11_cohort_treatment_naive_monocytes.ipynb) - discovery cohort monocyte analysis
14. [`14_starcat.ipynb`](14_starcat.ipynb) - starCAT T cell gene program projection into validation datasets
15. [`15_cohort_nonresponders.ipynb`](15_cohort_nonresponders.ipynb) - non-responder analysis
16. [`16_csf_validation.ipynb`](16_csf_validation.ipynb) - CSF/blood validation analysis
17. [`17_brain_tissue_validation.ipynb`](17_brain_tissue_validation.ipynb) - brain tissue validation analysis
18. [`18_natalizumab_validation.ipynb`](18_natalizumab_validation.ipynb) - natalizumab validation analysis