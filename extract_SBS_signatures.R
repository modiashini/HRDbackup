library(maftools)
library(tidyverse)
library(sigminer)
library(CompressiveNMF)

# Matching function
match_to_RefSigs <- function(sigs, ref){
  all <- sigminer::cosine(sigs, ref)
  ids <- cbind(1:nrow(all), apply(all, 1, which.max))
  data.frame("sigs" = colnames(sigs), 
             "best" = colnames(ref)[ids[, 2]], 
             "cosine" = round(all[ids], 3))
}

# Load MAF
df_maf <- readRDS("/Users/ashinimodi/Documents/Breast_AdenoCa_ICGC_maf.rds.gzip")
maf <- read.maf(df_maf)

# Generate SBS matrix
MutMatrices <- sig_tally(maf, mode = "ALL")
X <- t(MutMatrices$SBS_96)

# Run CompressiveNMF MAP
SigPrior <- 2000 * CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37 + 1
set.seed(42)
res <- CompressiveNMF::CompressiveNMF_map(X, S = SigPrior, K = 10, a = 1, alpha = 1)

# Match to COSMIC
match_table <- match_to_RefSigs(res$Signatures, CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37)

# Save
write.csv(match_table, file = "/Users/ashinimodi/Documents/SBS_cosine_matches.csv", row.names = FALSE)
write.csv(res$Signatures, file = "/Users/ashinimodi/Documents/SBS_signatures.csv")
write.csv(res$Theta, file = "/Users/ashinimodi/Documents/SBS_exposures.csv")
