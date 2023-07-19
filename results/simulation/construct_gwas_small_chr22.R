########################################################### pop1 ###########################################################
library(mvtnorm)
library(susieR)
library(data.table)
library(RhpcBLASctl)
set.seed(10)
blas_set_num_threads(10)
source("/import/home/share/mingxuan/fine_mapping/utils.R")
source("/import/home/share/mingxuan/fine_mapping/XMAP.R")
bim <- fread("/import/home/share/mingxuan/wegene_qc2/height_ukb_qc2_chr22.bim")
bp_chr <- bim$V4

nrep <- 50

# n1_all <- c(5000, 10000, 15000, 20000)
n1 <- 5000
K_true <- 3

h1_omega <- 0.005
h2_omega <- 0.005
rho_omega <- 0.8

enrich <- 50
rho_sigma <- 0

loci <- fread(paste0("/import/home/share/mingxuan/fine_mapping/loci_chr22_1mb.txt"))
i_loci <- 29
idx_include <- (501:1000) + loci$idx_start[i_loci] - 1
p <- length(idx_include)
p_chr <- nrow(bim)

# per-snp variance for bg
omega1 <- h1_omega / p
omega2 <- h2_omega / p

# per-snp variance for causal
Sigma1 <- omega1 * enrich
Sigma2 <- omega2 * enrich

# residual variance
Sig_E <- matrix(c(1 - omega1 * p - Sigma1 * K_true, 0, 0, 1 - omega2 * p - Sigma2 * K_true), 2, 2)



Beta1 <- data.matrix(fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_loci_",i_loci,"_small_p",p,"_Beta1_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                  "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".txt")))


Phi1 <- data.matrix(fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_loci_",i_loci,"_small_p",p,"_Phi1_chr22_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                 "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".txt")))


X1 <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_X1_wg_in_sample_ref_n", n1, "_chr22_loci", i_loci, ".txt"))
Xs1 <- scale(data.matrix(X1[, -1])[, idx_include])

fam <- fread("/import/home/share/mingxuan/wegene_qc2/height_merge_qc2_chr22.fam")
X1_22 <- read_data("/import/home/share/mingxuan/wegene_qc2/height_merge_qc2_chr22")

inds <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_indiv_wg_in_sample_ref_n", n1, ".txt"))
inds <- fam$V1[fam$V1 %in% inds$x]
X1_22 <- X1_22$X[match(inds, fam$V1),]
Xs1_22 <- scale(X1_22)

for (i in 1:nrep) {

  y1 <- Xs1 %*% Beta1[, i] +
    Xs1_22 %*% Phi1[, i] +
    rnorm(n1, 0, sqrt(1 - omega1 * p_chr - Sigma1 * K_true))

  rownames(y1) <- X1$V1

  sumstats1 <- univariate_regression(X1_22, y1)
  z_scores1 <- sumstats1$betahat / sumstats1$sebetahat
  # susie_plot(z_scores1,"z",b=Beta1[,i])

  names(z_scores1) <- bim$V2

  write.table(y1, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_y1_wholeCHR_loci_",i_loci,"_small_p",p,"_chr22_n", n1, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"),
              col.names = F, row.names = T, quote = F, sep = "\t")
  write.table(z_scores1, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_zs1_wholeCHR_loci_",i_loci,"_small_p",p,"_chr22_n", n1, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                       "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"),
              col.names = F, row.names = T, quote = F, sep = "\t")


  cat(i, "-th rep finished.\n")
}


########################################################### pop2 ###########################################################
library(mvtnorm)
library(susieR)
library(data.table)
library(RhpcBLASctl)
set.seed(1)
blas_set_num_threads(30)
source("/import/home/share/mingxuan/fine_mapping/utils.R")
source("/import/home/share/mingxuan/fine_mapping/XMAP.R")
bim <- fread("/import/home/share/mingxuan/wegene_qc2/height_ukb_qc2_chr22.bim")
bp_chr <- bim$V4

nrep <- 50

# n2_all <- c(10000, 15000, 20000)
n2 <- 20000
K_true <- 1

h1_omega <- 0.005
h2_omega <- 0.005
rho_omega <- 0.8

enrich <- 50
rho_sigma <- 0

loci <- fread(paste0("/import/home/share/mingxuan/fine_mapping/loci_chr22_1mb.txt"))
i_loci <- 29
idx_include <- (501:1000) + loci$idx_start[i_loci] - 1
p <- length(idx_include)
p_chr <- nrow(bim)


# per-snp variance for bg
omega1 <- h1_omega / p
omega2 <- h2_omega / p

# per-snp variance for causal
Sigma1 <- omega1 * enrich
Sigma2 <- omega2 * enrich

# residual variance
Sig_E <- matrix(c(1 - omega1 * p - Sigma1 * K_true, 0, 0, 1 - omega2 * p - Sigma2 * K_true), 2, 2)

Beta2 <- data.matrix(fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_loci_",i_loci,"_small_p",p,"_Beta2_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                  "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".txt")))

Phi2 <- data.matrix(fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_loci_",i_loci,"_small_p",p,"_Phi2_chr22_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                 "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".txt")))

X2 <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_X2_ukb_in_sample_ref_n", n2, "_chr22_loci", i_loci, ".txt"))
Xs2 <- scale(data.matrix(X2[, -1])[,idx_include])

fam <- fread("/import/home/share/mingxuan/wegene_qc2/height_ukb_qc2_chr22.fam")
X2_22 <- read_data("/import/home/share/mingxuan/wegene_qc2/height_ukb_qc2_chr22")

inds <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_indiv_ukb_in_sample_ref_n", n2, ".txt"))
inds <- fam$V1[fam$V1 %in% inds$x]
X2_22 <- X2_22$X[match(inds, fam$V1),]
Xs2_22 <- scale(X2_22)

for (i in 1:nrep) {

  y2 <- Xs2 %*% Beta2[, i] +
    Xs2_22 %*% Phi2[, i] +
    rnorm(n2, 0, sqrt(1 - omega2 * p_chr - Sigma2 * K_true))

  rownames(y2) <- X2$V1

  sumstats2 <- univariate_regression(X2_22, y2)
  z_scores2 <- sumstats2$betahat / sumstats2$sebetahat

  names(z_scores2) <- bim$V2

  write.table(y2, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_y2_wholeCHR_loci_",i_loci,"_small_p",p,"_chr22_n", n2, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"),
              col.names = F, row.names = T, quote = F, sep = "\t")
  write.table(z_scores2, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_zs2_wholeCHR_loci_",i_loci,"_small_p",p,"_chr22_n", n2, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                       "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"),
              col.names = F, row.names = T, quote = F, sep = "\t")
  cat(i, "-th rep finished.\n")
}



