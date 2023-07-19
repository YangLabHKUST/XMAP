library(VCM)
library(XPASS)
library(mvtnorm)
library(susieR)
library(data.table)
library(RhpcBLASctl)
set.seed(10)
blas_set_num_threads(10)
source("/import/home/share/mingxuan/fine_mapping/utils.R")
source("/import/home/share/mingxuan/fine_mapping/XMAP.R")
source("/home/share/mingxuan/taam/code/sldxr.R")
bim <- fread("/home/share/mingxuan/wegene_qc2/height_ukb_qc2_chr22.bim")
bp_chr <- bim$V4

nrep <- 50

# n1_all <- c(5000, 10000, 15000, 20000)
n1 <- 20000
n2 <- 20000
K_true <- 3


rho_sigma <- 0

pc_var1 <- 0.05
pc_var2 <- 0.2

i_loci <- 3
loci <- fread(paste0("/home/share/mingxuan/fine_mapping/loci_chr22_1mb.txt"))
idx_include <- loci$idx_start[i_loci]:loci$idx_end[i_loci] #(501:1000) + loci$idx_start[i_loci] - 1
p <- length(idx_include)

# per-snp variance for causal
Sigma1 <- 0.0005
Sigma2 <- 0.0005
Sigma12 <- rho_sigma * sqrt(Sigma1 * Sigma2)
Sigma <- matrix(c(Sigma1, Sigma12, Sigma12, Sigma2), 2, 2)

idx_candidate <- c(661, 800, 1300)



Beta1 <- Beta2 <- matrix(NA, p, nrep)
for (i in 1:nrep) {


  idx_causal <- sample(idx_candidate, K_true)
  idx_causal <- idx_causal[order(idx_causal)]


  beta <- matrix(0, p, 2)
  beta[idx_causal,] <- rmvnorm(K_true, rep(0, 2), Sigma)

  Beta1[, i] <- beta[, 1] / sqrt(var(beta[, 1])) * sqrt(Sigma1 / p)
  Beta2[, i] <- beta[, 2] / sqrt(var(beta[, 2])) * sqrt(Sigma2 / p)


  cat(i, "-th rep finished.\n")
}


write.table(Beta1, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_loci_",i_loci,"_p",p,"_Beta1_Ktrue", K_true,
                                 "_no_background_Sigma", Sigma1, "_", Sigma2, "_rhoSigma", rho_sigma, "_nrep", nrep, ".txt"),
            col.names = T, row.names = F, quote = F, sep = "\t")
write.table(Beta2, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_loci_",i_loci,"_p",p,"_Beta2_Ktrue", K_true,
                                 "_no_background_Sigma", Sigma1, "_", Sigma2, "_rhoSigma", rho_sigma, "_nrep", nrep, ".txt"),
            col.names = T, row.names = F, quote = F, sep = "\t")


