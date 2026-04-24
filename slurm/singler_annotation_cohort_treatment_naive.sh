#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=singler_annotation
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=1:00:00

BASEDIR="/project/fuggerlab/rfarooq/project/ocrelizumab_multiome" # change as needed

mkdir -p "$BASEDIR/results/tables/cite_seq/cohort_treatment_naive/cluster_annotate"

"$BASEDIR/scripts/singler_annotation.R" \
    --query "$BASEDIR/data/processed/cite_seq/cohort_treatment_naive/multi/peak-method:None_b:False_d:False/n-features:5000_normalisation:LogNormalize_clr:original_M:True_C:True/normalised.qs" \
    --reference MonacoImmuneData \
    --vars 'label.main;label.fine' \
    --output "$BASEDIR/results/tables/cite_seq/cohort_treatment_naive/cluster_annotate/singler_annotation.tsv" \
    --threads 4 \
    --log "$BASEDIR/results/tables/cite_seq/cohort_treatment_naive/cluster_annotate/singler_annotation.log"