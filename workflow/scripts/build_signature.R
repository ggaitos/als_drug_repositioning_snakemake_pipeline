suppressMessages({
  library(data.table); library(tidyverse); library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
spx_file <- args[1]   # spinal cord S-PrediXcan CSV
ann_file <- args[2]   # gencode table
out_dir  <- args[3]   # output directory

data <- fread(spx_file, header = TRUE)

# FDR adjust
Padj.FDR <- p.adjust(data$pvalue, method = "fdr")
print <- cbind(data, Padj.FDR)
write.csv(print, file.path(out_dir, "VR.ALS.quick.harmo_Brain_Spinal_cord_cervical_c-1.FDR.csv"))

# Filter to significant, split UP/DN
print <- filter(print, Padj.FDR < 0.05)
print$Signature[print$zscore > 0] <- "UP"
print$Signature[print$zscore < 0] <- "DN"

UP <- filter(print, zscore > 0); cat("UP genes:", nrow(UP), "\n")
DN <- filter(print, zscore < 0); cat("DN genes:", nrow(DN), "\n")

# The signature file that feeds SignatureSearch
sig <- dplyr::select(print, gene, gene_name, zscore, pvalue, Padj.FDR, Signature)
write.table(sig, file.path(out_dir, "ALS.SpinalCord.Signature"),
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# Annotated version (gene -> chr/pos) for the Miami plot
ann <- fread(ann_file, header = TRUE)
ann <- dplyr::select(ann, gene = Geneid, gene_name = GeneSymbol,
                     Chromosome, Start, End, Class)
ann$chr <- as.numeric(gsub("chr", "", ann$Chromosome))

fdr <- read.csv(file.path(out_dir, "VR.ALS.quick.harmo_Brain_Spinal_cord_cervical_c-1.FDR.csv"), header = TRUE)
merged <- merge(ann, fdr, by = "gene_name")
write.table(merged, file.path(out_dir, "ALL.ALS.SpinalCord.Signature.Ann.ENCODEv41.hg38.tab"),
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

cat("done\n")