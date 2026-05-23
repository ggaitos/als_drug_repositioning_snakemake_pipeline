#!/bin/bash
#SBATCH --account=PAS2598
#SBATCH --job-name=als_impute
#SBATCH --array=0-219
#SBATCH --time=2:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=results/ALS_VanRheenen2021/summary_imputation/logs/impute_%a.log

cd /fs/scratch/PDE0075/als-drug-repositioning

# 220 tasks = 22 chromosomes x 10 sub-batches
# task id -> chromosome (1..22) and sub-batch (0..9)
CHR=$(( SLURM_ARRAY_TASK_ID / 10 + 1 ))
SB=$(( SLURM_ARRAY_TASK_ID % 10 ))

envs/imlabtools/bin/python workflow/external/summary-gwas-imputation/src/gwas_summary_imputation.py \
  -by_region_file resources/data/eur_ld.bed.gz \
  -gwas_file results/ALS_VanRheenen2021/harmonized_gwas/ALS.VanRheenen2021_quick.Formatted.harmonized.tabul.txt.gz \
  -parquet_genotype resources/data/reference_panel_1000G/chr${CHR}.variants.parquet \
  -parquet_genotype_metadata resources/data/reference_panel_1000G/variant_metadata.parquet \
  -window 100000 \
  -parsimony 7 \
  -chromosome ${CHR} \
  -regularization 0.1 \
  -frequency_filter 0.01 \
  -sub_batches 10 \
  -sub_batch ${SB} \
  --standardise_dosages \
  -output results/ALS_VanRheenen2021/summary_imputation/VR.ALS.quick.harmo.YESfreq_chr${CHR}_sb${SB}.txt.gz