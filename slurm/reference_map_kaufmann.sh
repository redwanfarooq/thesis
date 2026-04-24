#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=reference_map
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --time=12:00:00

BASEDIR="/project/fuggerlab/rfarooq/project/ocrelizumab_multiome" # change as needed

mkdir -p "$BASEDIR/results/tables/external/kaufmann/cluster_annotate"

"$BASEDIR/scripts/reference_map.R" \
    --query "$BASEDIR/data/processed/external/kaufmann/annotated/kaufmann_full.qs" \
    --reference "$BASEDIR/data/ref/datasets/pbmc_ocrelizumab_cohort_treatment_naive.qs" \
    --vars 'cluster_main;cluster_coarse;cluster_fine' \
    --output "$BASEDIR/results/tables/external/kaufmann/cluster_annotate/reference_map.tsv" \
    --threads 4 \
    --log "$BASEDIR/results/tables/external/kaufmann/cluster_annotate/reference_map.log" \
    --exp HASH \
    --query-assay RNA \
    --reference-assay RNA \
    --reference-reduction spca \
    --normalisation-method LogNormalize