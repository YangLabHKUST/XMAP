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
blas_set_num_threads(30)
source("/import/home/share/mingxuan/fine_mapping/utils.R")
source("/import/home/share/mingxuan/fine_mapping/XMAP.R")
source("/import/home/share/mingxuan/taam/code/sldxr.R")
bim <- fread("/import/home/share/mingxuan/wegene_qc2/height_ukb_qc2_chr22.bim")
bp_chr <- bim$V4

nrep <- 50

# n1_all <- c(5000, 10000, 15000, 20000)
n1 <- 15000#as.numeric(args[1]) #5000
K_true <- 3#as.numeric(args[2]) #1

h1_omega <- 0.005
h2_omega <- 0.005
rho_omega <- 0.8

enrich <- 50
rho_sigma <- 0

pc_var <- 0.05#as.numeric(args[3]) #0.1

i_loci <- 29
loci <- fread(paste0("/import/home/share/mingxuan/fine_mapping/loci_chr22_1mb.txt"))
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

fam <- fread("/import/home/share/mingxuan/wegene_qc2/height_merge_qc2_chr22.fam")
X1_22 <- read_data("/import/home/share/mingxuan/wegene_qc2/height_merge_qc2_chr22")

inds <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_indiv_wg_in_sample_ref_n", n1, ".txt"))
pcs <- fread("/import/home/share/xiaojs/height/genodata/height_weight_pc_ances_covar_20190924.txt")
inds <- fam$V1[fam$V1 %in% inds$x]
pc1 <- pcs$pc1[match(inds, pcs$barcode)]
X1_22 <- X1_22$X[match(inds, fam$V1),]

for (i in 1:nrep) {


  y1 <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_y1_inflate_pc1_", pc_var, "_loci_", i_loci, "_small_p", p, "_chr22_n", n1, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                     "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"))
  z_scores1 <- rep(NA, ncol(X1_22))
  for (j in 1:ncol(X1_22)) {
    z_scores1[j] <- summary(lm(y ~ ., data.frame(y = y1$V2, x = X1_22[, j], pc = pc1)))$coefficients[2, 3]
    if (j %% 1000 == 0) cat(i, "-th rep, ", j, "/", ncol(X1_22), " SNP finished.\n")
  }
  names(z_scores1) <- bim$V2
  write.table(z_scores1, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_zs1_inflate_pc1_", pc_var, "_corrected_loci_", i_loci, "_small_p", p, "_chr22_n", n1, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
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
cores = 30
cl <- makeCluster(cores) #not to overload your computer
registerDoParallel(cl)
set.seed(1)
blas_set_num_threads(cores)
source("/import/home/share/mingxuan/fine_mapping/utils.R")
source("/import/home/share/mingxuan/fine_mapping/XMAP.R")
bim <- fread("/import/home/share/mingxuan/wegene_qc2/height_ukb_qc2_chr22.bim")
bp_chr <- bim$V4

nrep <- 50

# n2_all <- c(10000, 15000, 20000)
n2 <- 20000 #as.numeric(args[1])
K_true <- 3 # as.numeric(args[2]) #3

h1_omega <- 0.005
h2_omega <- 0.005
rho_omega <- 0.8

enrich <- 50
rho_sigma <- 0

pc_var <- 0.2 # as.numeric(args[3])

loci <- fread(paste0("/import/home/share/mingxuan/fine_mapping/loci_chr22_1mb.txt"))
i_loci <- 29
idx_include <- (501:1000) + loci$idx_start[i_loci] - 1
p <- length(idx_include)


# per-snp variance for bg
omega1 <- h1_omega / p
omega2 <- h2_omega / p

# per-snp variance for causal
Sigma1 <- omega1 * enrich
Sigma2 <- omega2 * enrich

Beta2 <- data.matrix(fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_loci_", i_loci, "_small_p", p, "_Beta2_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                  "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".txt")))


X2 <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_X2_ukb_in_sample_ref_n", n2, "_chr22_loci", i_loci, ".txt"))
Xs2 <- scale(data.matrix(X2[, -1])[, idx_include])

fam <- fread("/import/home/share/mingxuan/wegene_qc2/height_ukb_qc2_chr22.fam")
X2_22 <- read_data("/import/home/share/mingxuan/wegene_qc2/height_ukb_qc2_chr22")

inds <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_indiv_ukb_in_sample_ref_n", n2, ".txt"))
inds <- fam$V1[fam$V1 %in% inds$x]
pcs <- fread("/import/home/share/xhu/database/UKB/ukb_sqc_v2_organized.txt")
pc1 <- pcs$pc1[match(inds, pcs$FID)]
X2_22 <- X2_22$X[match(inds, fam$V1),]

for (i in 1:nrep) {


  y2 <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_y2_inflate_pc1_", pc_var, "_loci_", i_loci, "_small_p", p, "_chr22_n", n2, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                     "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"))

  z_scores2 <- rep(NA, ncol(X2_22))
  for (j in 1:ncol(X2_22)) {
    z_scores2[j] <- summary(lm(y ~ ., data.frame(y = y2$V2, x = X2_22[, j], pc = pc1)))$coefficients[2, 3]
    if (j %% 1000 == 0) cat(i, "-th rep, ", j, "/", ncol(X2_22), " SNP finished.\n")
  }

  names(z_scores2) <- bim$V2
  write.table(z_scores2, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_zs2_inflate_pc1_", pc_var, "_corrected_loci_", i_loci, "_small_p", p, "_chr22_n", n2, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                       "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"),
              col.names = F, row.names = T, quote = F, sep = "\t")
  cat(i, "-th rep finished.\n")
}


