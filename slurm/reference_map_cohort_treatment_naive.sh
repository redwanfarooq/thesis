#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=reference_map
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --time=12:00:00

BASEDIR="/project/fuggerlab/rfarooq/project/ocrelizumab_multiome" # change as needed

mkdir -p "$BASEDIR/results/tables/cite_seq/cohort_treatment_naive/cluster_annotate"

"$BASEDIR/scripts/reference_map.R" \
    --query "$BASEDIR/data/processed/cite_seq/cohort_treatment_naive/multi/peak-method:None_b:False_d:False/n-features:5000_normalisation:LogNormalize_clr:seurat_M:True_C:True/normalised.qs" \
    --reference "$BASEDIR/data/ref/datasets/pbmc_multimodal_2023.rds" \
    --vars 'celltype.l1;celltype.l2;celltype.l3' \
    --output "$BASEDIR/results/tables/cite_seq/cohort_treatment_naive/cluster_annotate/reference_map.tsv" \
    --threads 4 \
    --log "$BASEDIR/results/tables/cite_seq/cohort_treatment_naive/cluster_annotate/reference_map.log" \
    --query-assay RNA \
    --reference-assay SCT \
    --normalisation-method SCT