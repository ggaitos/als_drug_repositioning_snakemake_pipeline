library(data.table)
library(tidyverse)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
harmo_file <- args[1]   # harmonized quick output
raw_file   <- args[2]   # original GWAS
out_file   <- args[3]   # output .txt (tab-delimited, pre-bgzip)

data <- fread(harmo_file, header = TRUE)
data$sample_size      <- 138086
data$n_cases          <- 27205
data$effect_size      <- data$beta
data$panel_variant_id <- data$snp
data$variant_id       <- data$gwas_snp
data <- data %>% separate(snp, c("chromosome", "position"), extra = "drop", fill = "right")

data2 <- fread(raw_file, header = TRUE)
data2$variant_id <- data2$rsid
data2$frequency  <- data2$effect_allele_frequency
data2 <- select(data2, variant_id, standard_error, frequency)

merged <- merge(data, data2, by = "variant_id")

dataprint <- select(merged, variant_id, panel_variant_id, chromosome, position,
                    effect_allele, non_effect_allele, frequency, pvalue, zscore,
                    effect_size, standard_error, sample_size, n_cases)
write.table(dataprint, out_file, quote = FALSE, row.names = FALSE, sep = "\t")