#!/bin/bash
#SBATCH --account=PDE0075
#SBATCH --job-name=als_spx
#SBATCH --array=0-13
#SBATCH --time=0:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=results/ALS_VanRheenen2021/spredixcan/logs/spx_%a.log

cd /fs/scratch/PDE0075/als-drug-repositioning

TISSUES=(Whole_Blood Brain_Amygdala Brain_Anterior_cingulate_cortex_BA24 \
  Brain_Caudate_basal_ganglia Brain_Cerebellar_Hemisphere Brain_Cerebellum \
  Brain_Cortex Brain_Frontal_Cortex_BA9 Brain_Hypothalamus \
  Brain_Nucleus_accumbens_basal_ganglia Brain_Putamen_basal_ganglia \
  Brain_Spinal_cord_cervical_c-1 Brain_Substantia_nigra)

T=${TISSUES[$SLURM_ARRAY_TASK_ID]}

envs/imlabtools/bin/python workflow/external/MetaXcan/software/SPrediXcan.py \
  --gwas_file results/ALS_VanRheenen2021/processed_summary_imputation/imputed_VR.ALS.quick.harmo.YESfreq.txt.gz \
  --snp_column panel_variant_id \
  --effect_allele_column effect_allele \
  --non_effect_allele_column non_effect_allele \
  --zscore_column zscore \
  --model_db_path resources/data/models/eqtl/mashr/mashr_${T}.db \
  --covariance resources/data/models/eqtl/mashr/mashr_${T}.txt.gz \
  --model_db_snp_key varID \
  --keep_non_rsid --additional_output --throw \
  --output_file results/ALS_VanRheenen2021/spredixcan/VR.ALS.quick.harmo.YESfreq_PM_${T}.csv