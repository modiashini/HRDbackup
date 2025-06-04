library(maftools)
library(tidyverse)
library(sigminer)
library(CompressiveNMF)

#----- Useful function to match to cosmic
match_to_RefSigs <- function(sigs, ref){
  all <- sigminer::cosine(sigs, ref)
  ids <- cbind(1:nrow(all), apply(all, 1, which.max))
  data.frame("sigs" = colnames(sigs), "best" = colnames(ref)[ids[, 2]], "cosine" = round(all[ids], 3))
}
#-----

df_maf <- readRDS("~/SigPoisProcess/Breast_AdenoCa_ICGC_maf.rds.gzip")
maf <- read.maf(df_maf)
MutMatrices <- sig_tally(maf, mode = "ALL")


# Indel signatures
Sid <- CompressiveNMF::COSMIC_v3.4_ID83_GRCh37
SigPrior <- 5000 * Sid + 1

# MAP solution
X <- t(MutMatrices$ID_83)
X <- X[rownames(Sid), ]  # Re-arrange the labels!
set.seed(42)
res <- CompressiveNMF::CompressiveNMF_map(X, S = SigPrior, K = 10, a = 1, alpha = 1)
CompressiveNMF::plot_ID_signature(res$Signatures)
match_to_RefSigs(res$Signatures, CompressiveNMF::COSMIC_v3.4_ID83_GRCh37)
CompressiveNMF::plot_weights(res$Theta)


# Bayesian solution
set.seed(42)
resBayes <- CompressiveNMF::CompressiveNMF(X,
                                      S = Sid,
                                      nsamples = 1000,
                                      burnin = 3500,
                                      betah_optimal = TRUE,
                                      K = 10,
                                      a = 0.5, alpha = 0.5,
                                      swap_prior = TRUE)

print(resBayes)
plot(resBayes)
CompressiveNMF::plot_weights(Wmat = resBayes$Weights, nclust = 10)
match_to_RefSigs(resBayes$Signatures, Sid)
plot(resBayes$Signatures %*% resBayes$Weights, X)



