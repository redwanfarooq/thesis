#!/bin/bash
#SBATCH --partition=long
#SBATCH --job-name=sctriangulate
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=600G
#SBATCH --time=3-00:00:00

BASEDIR="/project/fuggerlab/rfarooq/project/ocrelizumab_multiome" # change as needed

mkdir -p "$BASEDIR/results/tables/cite_seq/cohort_treatment_naive/cluster_annotate/sctriangulate"

"$BASEDIR/scripts/run_sctriangulate.py" \
    --input "$BASEDIR/results/tables/cite_seq/cohort_treatment_naive/cluster_annotate/clustered_annotated.h5mu" \
    --outdir "$BASEDIR/results/tables/cite_seq/cohort_treatment_naive/cluster_annotate/sctriangulate" \
    --query 'RNA:leiden_2;RNA:leiden_5;RNA:leiden_8;ADT:leiden_2;ADT:leiden_5;celltypist_L3;azimuth_celltype.l3;singler_label.fine' \
    --threads 16 \
    --log "$BASEDIR/results/tables/cite_seq/cohort_treatment_naive/cluster_annotate/sctriangulate/run_sctriangulate.log" \
    --counts-layer counts \
    --compute-shapley-parallel
