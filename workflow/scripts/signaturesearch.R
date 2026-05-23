suppressMessages({
  library(data.table); library(tidyverse); library(dplyr)
  library(org.Hs.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)
spx_file <- args[1]   # spinal cord S-PrediXcan CSV
out_dir  <- args[2]
eh_cache <- args[3]   # ExperimentHub cache dir

Sys.setenv(EXPERIMENT_HUB_CACHE = eh_cache)

# --- FDR, filter, UP/DN (same as cell 5) ---
data <- fread(spx_file, header = TRUE)
Padj.FDR <- p.adjust(data$pvalue, method = "fdr")
print <- cbind(data, Padj.FDR)
print <- filter(print, Padj.FDR < 0.05)

# Ensembl (strip version) -> Entrez
print <- print %>% separate(gene, c("ENS"), extra = "drop", fill = "right")
con <- mapIds(org.Hs.eg.db, as.character(print$ENS), "ENTREZID", "ENSEMBL")
print$geneID <- con[match(print$ENS, names(con))]
print <- filter(print, !is.na(geneID))

upset   <- as.character(filter(print, zscore > 0)$geneID)
downset <- as.character(filter(print, zscore < 0)$geneID)
cat("upset:", length(upset), " downset:", length(downset), "\n")

# --- Drug search (cell 5 core) ---
suppressMessages({
  library(signatureSearch); library(ExperimentHub); library(rhdf5)
  library(mygene); library(Hmisc)
})
eh <- ExperimentHub()
lincs <- eh[["EH3226"]]

# LINCS method on LINCS db -- this is the one that feeds the ranking
qsig_lincs <- qSig(query = list(upset = upset, downset = downset),
                   gess_method = "LINCS", refdb = lincs)
lincsm <- gess_lincs(qsig_lincs, sortby = "NCS", tau = TRUE, workers = 1)
write.csv(as.data.frame(lincsm@result),
          file.path(out_dir, "VR.ALS.lincsmethod.lincsDS.SpinalCordc1FDR.csv"))

# CMAP method on LINCS db (comparator)
qsig_cmap <- qSig(query = list(upset = upset, downset = downset),
                  gess_method = "CMAP", refdb = lincs)
cmapm <- gess_cmap(qSig = qsig_cmap, chunk_size = 5000, workers = 1)
write.csv(as.data.frame(cmapm@result),
          file.path(out_dir, "VR.ALS.cmapmethod.lincsDS.SpinalCordc1FDR.csv"))

cat("done\n")