library(mvtnorm)
library(susieR)
library(data.table)
library(RhpcBLASctl)
set.seed(45)#set.seed(5)
blas_set_num_threads(20)
source("/import/home/share/mingxuan/fine_mapping/utils.R")
source("/import/home/share/mingxuan/fine_mapping/XMAP.R")

bim_chr <- fread("/import/home/share/mingxuan/wegene_qc2/height_ukb_qc2_chr22.bim")
bp_chr <- bim_chr

nrep <- 50


K_true <- 3

h1_omega <- 0.005
h2_omega <- 0.005
rho_omega <- 0.8

enrich <- 50
rho_sigma <- 0


loci <- fread(paste0("/import/home/share/mingxuan/fine_mapping/loci_chr22_1mb.txt"))


# candidate SNPs selected by #(r_ukb>0.9)>=3 and #(r_wg>0.6)=0 in loci 29
i_loci <- 29
idx_include <- (501:1000) + loci$idx_start[i_loci] - 1
p <- length(idx_include)
p_chr <- nrow(bim_chr)

idx_candidate <- c(709,780,886) - 501 + 1

# per-snp variance for bg
omega1 <- h1_omega / p
omega2 <- h2_omega / p
omega12 <- rho_omega * sqrt(omega1 * omega2)

Omega <- matrix(c(omega1, omega12, omega12, omega2), 2, 2)


# per-snp variance for causal
Sigma1 <- omega1 * enrich
Sigma2 <- omega2 * enrich
Sigma12 <- rho_sigma * sqrt(Sigma1 * Sigma2)
Sigma <- matrix(c(Sigma1, Sigma12, Sigma12, Sigma2), 2, 2)

Phi1 <- Phi2 <- matrix(NA, p_chr, nrep)


Beta1 <- Beta2 <- matrix(NA, p, nrep)

for (i in 1:nrep) {


  phi <- rmvnorm(p_chr, rep(0, 2), Omega)

  idx_causal <- sample(idx_candidate, K_true)
  idx_causal <- idx_causal[order(idx_causal)]


  beta <- matrix(0, p, 2)
  beta[idx_causal,] <- rmvnorm(K_true, rep(0, 2), Sigma)

  Beta1[, i] <- beta[, 1]/sqrt(var(beta[,1]))*sqrt(Sigma[1,1]/p*K_true)
  Beta2[, i] <- beta[, 2]/sqrt(var(beta[,2]))*sqrt(Sigma[2,2]/p*K_true)


  Phi1[, i] <- phi[, 1]
  Phi2[, i] <- phi[, 2]


  cat(i, "-th rep finished.\n")
}
write.table(Beta1, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_loci_",i_loci,"_small_p",p,"_Beta1_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                 "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".txt"),
            col.names = T, row.names = F, quote = F, sep = "\t")
write.table(Beta2, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_loci_",i_loci,"_small_p",p,"_Beta2_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                 "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".txt"),
            col.names = T, row.names = F, quote = F, sep = "\t")


write.table(Phi1, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_loci_",i_loci,"_small_p",p,"_Phi1_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".txt"),
            col.names = T, row.names = F, quote = F, sep = "\t")

write.table(Phi2, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_loci_",i_loci,"_small_p",p,"_Phi2_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".txt"),
            col.names = T, row.names = F, quote = F, sep = "\t")



