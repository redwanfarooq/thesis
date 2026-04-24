#!/bin/bash
#SBATCH --partition=long
#SBATCH --job-name=meld_simulations
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=150G
#SBATCH --time=3-00:00:00

BASEDIR="/project/fuggerlab/rfarooq/project/ocrelizumab_multiome" # change as needed

"$BASEDIR/scripts/meld_simulations.py" \
    --pca "$BASEDIR/data/processed/cite_seq/cohort_treatment_naive/embeddings/integrated_rna.tsv.gz" \
    --output "$BASEDIR/results/tables/cite_seq/cohort_treatment_naive/meld_simulations" \
    --log "$BASEDIR/results/tables/cite_seq/cohort_treatment_naive/meld_simulations/meld_simulations.log" \
    --metadata "$BASEDIR/results/tables/cite_seq/cohort_treatment_naive/meld_simulations/metadata.tsv" \
    --replicate donor \
    --condition condition \
    --reference Baseline \
    --k 25 \
    --beta 90 \
    --n-iter 10000 \
    -t 25

# k and beta parameters chosen based on hyperparameter grid search