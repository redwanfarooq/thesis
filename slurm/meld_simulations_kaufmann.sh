#!/bin/bash
#SBATCH --partition=long
#SBATCH --job-name=meld_simulations
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=150G
#SBATCH --time=3-00:00:00

BASEDIR="/project/fuggerlab/rfarooq/project/ocrelizumab_multiome" # change as needed

"$BASEDIR/scripts/meld_simulations.py" \
    --pca "$BASEDIR/data/processed/external/kaufmann/embeddings/integrated_rna.tsv.gz" \
    --output "$BASEDIR/results/tables/external/kaufmann/meld_simulations" \
    --log "$BASEDIR/results/tables/external/kaufmann/meld_simulations/meld_simulations.log" \
    --metadata "$BASEDIR/results/tables/external/kaufmann/meld_simulations/metadata.tsv" \
    --replicate donor \
    --condition condition \
    --reference Pre \
    --k 20 \
    --beta 70 \
    --n-iter 10000 \
    -t 25

# k and beta parameters chosen based on hyperparameter grid search