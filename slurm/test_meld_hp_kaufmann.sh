#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=test_meld_hp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=250G
#SBATCH --time=1-00:00:00

BASEDIR="/project/fuggerlab/rfarooq/project/ocrelizumab_multiome" # change as needed

"$BASEDIR/scripts/test_meld_hp.py" \
    --pca "$BASEDIR/data/processed/external/kaufmann/embeddings/integrated_rna.tsv.gz" \
    --output "$BASEDIR/results/tables/external/kaufmann/test_meld_hp.tsv" \
    --log "$BASEDIR/results/tables/external/kaufmann/test_meld_hp.log" \
    --k-lower 5 \
    --k-upper 51 \
    --k-step 5 \
    --b-lower 10 \
    --b-upper 201 \
    --b-step 10 \
    -t 25