library(data.table); library(dplyr)
args <- commandArgs(trailingOnly=TRUE)
lincs_file <- args[1]; approved_file <- args[2]; out_file <- args[3]

approved <- trimws(tolower(readLines(approved_file)))
lincs <- fread(lincs_file); lincs$pert <- tolower(lincs$pert)
fne <- function(x){ x<-x[!is.na(x)&x!=""&x!="NA"]; if(length(x)) x[1] else NA }

df <- lincs %>% filter(pert %in% approved) %>%
  group_by(pert) %>%
  summarise(n_cells=n(), NCS=mean(NCS,na.rm=TRUE), Tau=mean(Tau,na.rm=TRUE),
            WTCS=mean(WTCS,na.rm=TRUE), MOAss=fne(MOAss),
            t_gn_sym=fne(t_gn_sym), PCIDss=fne(PCIDss), .groups="drop") %>%
  as.data.frame()
df <- df[order(df$NCS),]
fwrite(df, out_file)
cat("approved drugs matched:", nrow(df), "\n")
