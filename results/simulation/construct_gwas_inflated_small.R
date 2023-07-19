#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
########################################################### pop1 ###########################################################
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
n1 <- as.numeric(args[1])#5000
K_true <- as.numeric(args[2])#1

h1_omega <- 0.005
h2_omega <- 0.005
rho_omega <- 0.8

enrich <- 50
rho_sigma <- 0

pc_var <- as.numeric(args[3])#0.1


i_loci <- 29
loci <- fread(paste0("/home/share/mingxuan/fine_mapping/loci_chr22_1mb.txt"))
idx_include <- (501:1000) + loci$idx_start[i_loci] - 1
p <- length(idx_include)

# per-snp variance for bg
omega1 <- h1_omega / p
omega2 <- h2_omega / p

# per-snp variance for causal
Sigma1 <- omega1 * enrich
Sigma2 <- omega2 * enrich




Beta1 <- data.matrix(fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_loci_", i_loci, "_small_p", p, "_Beta1_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                  "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".txt")))


X1 <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_X1_wg_in_sample_ref_n", n1, "_chr22_loci", i_loci, ".txt"))
Xs1 <- scale(data.matrix(X1[, -1])[, idx_include])

fam <- fread("/home/share/mingxuan/wegene_qc2/height_merge_qc2_chr22.fam")
X1_22 <- read_data("/home/share/mingxuan/wegene_qc2/height_merge_qc2_chr22")

inds <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_indiv_wg_in_sample_ref_n", n1, ".txt"))
pcs <- fread("/home/share/xiaojs/height/genodata/height_weight_pc_ances_covar_20190924.txt")
inds <- fam$V1[fam$V1 %in% inds$x]
pc1 <- pcs$pc1[match(inds, pcs$barcode)]
X1_22 <- X1_22$X[match(inds, fam$V1),]

for (i in 1:nrep) {


  inflate_term <- scale(pc1) * sqrt(pc_var)
  y1 <- inflate_term +
    Xs1 %*% Beta1[, i] +
    rnorm(n1, 0, sqrt(1 - Sigma1 * K_true - pc_var))

  rownames(y1) <- X1$V1

  sumstats1 <- univariate_regression(X1_22, y1)
  z_scores1 <- sumstats1$betahat / sumstats1$sebetahat


  names(z_scores1) <- bim$V2

  write.table(y1, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_y1_inflate_pc1_",pc_var,"_loci_", i_loci, "_small_p", p, "_chr22_n", n1, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"),
              col.names = F, row.names = T, quote = F, sep = "\t")
  write.table(z_scores1, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_zs1_inflate_pc1_",pc_var,"_loci_", i_loci, "_small_p", p, "_chr22_n", n1, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                       "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"),
              col.names = F, row.names = T, quote = F, sep = "\t")


  cat(i, "-th rep finished.\n")
}


########################################################### pop2 ###########################################################
library(VCM)
library(XPASS)
library(mvtnorm)
library(susieR)
library(data.table)
library(RhpcBLASctl)
set.seed(1)
blas_set_num_threads(10)
source("/import/home/share/mingxuan/fine_mapping/utils.R")
source("/import/home/share/mingxuan/fine_mapping/XMAP.R")
bim <- fread("/home/share/mingxuan/wegene_qc2/height_ukb_qc2_chr22.bim")
bp_chr <- bim$V4

nrep <- 50

# n2_all <- c(10000, 15000, 20000)
n2 <- as.numeric(args[1])#20000
K_true <- as.numeric(args[2])#3

h1_omega <- 0.005
h2_omega <- 0.005
rho_omega <- 0.8

enrich <- 50
rho_sigma <- 0

pc_var <- as.numeric(args[3])#0.2

loci <- fread(paste0("/home/share/mingxuan/fine_mapping/loci_chr22_1mb.txt"))
i_loci <- 29
idx_include <- (501:1000) + loci$idx_start[i_loci] - 1
p <- length(idx_include)


# per-snp variance for bg
omega1 <- h1_omega / p
omega2 <- h2_omega / p

# per-snp variance for causal
Sigma1 <- omega1 * enrich
Sigma2 <- omega2 * enrich

# residual variance

Beta2 <- data.matrix(fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_loci_", i_loci, "_small_p", p, "_Beta2_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                  "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".txt")))


X2 <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_X2_ukb_in_sample_ref_n", n2, "_chr22_loci", i_loci, ".txt"))
Xs2 <- scale(data.matrix(X2[, -1])[, idx_include])

fam <- fread("/home/share/mingxuan/wegene_qc2/height_ukb_qc2_chr22.fam")
X2_22 <- read_data("/home/share/mingxuan/wegene_qc2/height_ukb_qc2_chr22")

inds <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_indiv_ukb_in_sample_ref_n", n2, ".txt"))
inds <- fam$V1[fam$V1 %in% inds$x]
pcs <- fread("/home/share/xhu/database/UKB/ukb_sqc_v2_organized.txt")
pc1 <- pcs$pc1[match(inds, pcs$FID)]
X2_22 <- X2_22$X[match(inds, fam$V1),]

for (i in 1:nrep) {
  inflate_term <- scale(pc1) * sqrt(pc_var)
  y2 <- inflate_term +
    Xs2 %*% Beta2[, i] +
    rnorm(n2, 0, sqrt(1 - Sigma2 * K_true - pc_var))

  rownames(y2) <- X2$V1

  sumstats2 <- univariate_regression(X2_22, y2)
  z_scores2 <- sumstats2$betahat / sumstats2$sebetahat
  susie_plot(z_scores2[loci$idx_left[i_loci]:loci$idx_right[i_loci]][idx_include],"z",b=Beta2[,i])

  names(z_scores2) <- bim$V2

  write.table(y2, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_y2_inflate_pc1_",pc_var,"_loci_", i_loci, "_small_p", p, "_chr22_n", n2, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"),
              col.names = F, row.names = T, quote = F, sep = "\t")
  write.table(z_scores2, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_zs2_inflate_pc1_",pc_var,"_loci_", i_loci, "_small_p", p, "_chr22_n", n2, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                       "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"),
              col.names = F, row.names = T, quote = F, sep = "\t")
  cat(i, "-th rep finished.\n")
}



