#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=convert_seurat_scanpy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=1:00:00

BASEDIR="/project/fuggerlab/rfarooq/project/ocrelizumab_multiome" # change as needed

"$BASEDIR/scripts/convert_seurat_scanpy.R" \
    --input "$BASEDIR/data/processed/cite_seq/cohort_treatment_naive/multi/peak-method:None_b:False_d:False/n-features:5000_normalisation:LogNormalize_clr:seurat_M:True_C:True/batch:orig.ident_normalisation:LogNormalize_integration:harmony/integrated.qs" \
    --outdir "$BASEDIR/data/processed/cite_seq/cohort_treatment_naive/multi/peak-method:None_b:False_d:False/n-features:5000_normalisation:LogNormalize_clr:seurat_M:True_C:True/batch:orig.ident_normalisation:LogNormalize_integration:harmony" \
    --threads 4 \
    --log "$BASEDIR/data/processed/cite_seq/cohort_treatment_naive/multi/peak-method:None_b:False_d:False/n-features:5000_normalisation:LogNormalize_clr:seurat_M:True_C:True/batch:orig.ident_normalisation:LogNormalize_integration:harmony/convert_seurat_scanpy.log" \
    --assays 'RNA;ADT' \
    --obsm 'RNA=integrated.rna;ADT=integrated.adt'