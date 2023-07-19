#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
########################################################### pop1 ###########################################################
library(XPASS)
library(mvtnorm)
library(susieR)
library(data.table)
library(RhpcBLASctl)
set.seed(1)
blas_set_num_threads(20)
source("/import/home/share/mingxuan/fine_mapping/utils.R")
source("/import/home/share/mingxuan/fine_mapping/XMAP.R")
source("/home/share/mingxuan/taam/code/sldxr.R")
bim <- fread("/home/share/mingxuan/wegene_qc2/height_ukb_qc2_chr22.bim")
bp_chr <- bim$V4

nrep <- 50

# n1_all <- c(5000, 10000, 15000, 20000)
# n2_all <- c(10000, 15000, 20000)
n1 <- 15000#as.numeric(args[1]) #10000
n2 <- 20000
K_true <- 1#as.numeric(args[2]) #1
K <- 5

h1_omega <- 0.005
h2_omega <- 0.005
rho_omega <- 0.8

enrich <- 50
rho_sigma <- 0

pc_var1 <- 0.05
pc_var2 <- 0.2

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

R1 <- data.matrix(fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_R1_wg_in_sample_ref_n", n1, "_chr22_loci", i_loci, ".txt")))
R1 <- R1[idx_include, idx_include]
Beta1 <- data.matrix(fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_loci_", i_loci, "_small_p", p, "_Beta1_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                  "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".txt")))
R2 <- data.matrix(fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_R2_ukb_in_sample_ref_n", n2, "_chr22_loci", i_loci, ".txt")))
R2 <- R2[idx_include, idx_include]
Beta2 <- data.matrix(fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_loci_", i_loci, "_small_p", p, "_Beta2_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                  "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".txt")))

ldscore <- fread(paste0("/home/share/mingxuan/fine_mapping/LD_mat/LDscore_UKB_wegene_chr22.txt"))

# output
cs_power <- cs_coverage <- power <- precision <- matrix(NA, nrep, 3)
out_xmap_CI <- out_xmap_C <- list()
out_ldscore <- list()
for (i in 1:nrep) {
  idx_causal <- which(Beta1[, i] != 0)
  z_scores1 <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_zs1_inflate_pc1_", pc_var1, "_loci_", i_loci, "_small_p", p, "_chr22_n", n1, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                            "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"))
  z_scores2 <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_zs2_inflate_pc1_", pc_var2, "_loci_", i_loci, "_small_p", p, "_chr22_n", n2, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                            "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"))
  # susie_plot(z_scores1$V2,"z",b=Beta1[, i])
  # susie_plot(z_scores2$V2,"z",b=Beta2[, i])

  # ldsc estimate parameters
  idx1 <- which(z_scores1$V2^2 < 30 & z_scores2$V2^2 < 30)
  ld1 <- ldscore$EAS
  ld2 <- ldscore$EUR
  ld12 <- ldscore$EUR_EAS

  ld_w1 <- 1 / sapply(ld1, function(x) max(x, 1))
  ld_w2 <- 1 / sapply(ld2, function(x) max(x, 1))
  fit_step1 <- estimate_gc(data.frame(Z = z_scores1$V2[idx1], N = n1), data.frame(Z = z_scores2$V2[idx1], N = n2),
                           ld1[idx1], ld2[idx1], ld12[idx1],
                           reg_w1 = ld_w1[idx1], reg_w2 = ld_w2[idx1], reg_wx = sqrt(ld_w1[idx1] * ld_w2[idx1]),
                           constrain_intercept = F)


  fit_step2 <- estimate_gc(data.frame(Z = z_scores1$V2, N = n1), data.frame(Z = z_scores2$V2, N = n2),
                           ld1, ld2, ld12,
                           reg_w1 = ld_w1, reg_w2 = ld_w2, reg_wx = sqrt(ld_w1 * ld_w2),
                           constrain_intercept = T, fit_step1$tau1[1], fit_step1$tau2[1], fit_step1$theta[1])
  # c(fit_step2$tau1[2] * p, fit_step2$tau2[2] * p, fit_step2$theta[2] / sqrt(fit_step2$tau1[2] * fit_step2$tau2[2]))

  out_ldscore[[i]] <- fit_step2
  # threshold bg estimate
  omega1Hat <- fit_step2$tau1[2]
  omega2Hat <- fit_step2$tau2[2]
  rhohat <- fit_step2$theta[2] / sqrt(fit_step2$tau1[2] * fit_step2$tau2[2])
  if (omega1Hat < 1e-30) {
    omega1Hat <- 1e-30
    rhohat <- 0
  }
  if (omega2Hat < 1e-30) {
    omega2Hat <- 1e-30
    rhohat <- 0
  }
  rhohat <- ifelse(rhohat > 1, 0.99, rhohat)
  rhohat <- ifelse(rhohat < -1, -0.99, rhohat)
  OmegaHat0 <- OmegaHat <- matrix(c(omega1Hat, 0, 0, omega2Hat), 2, 2)
  OmegaHat[1, 2] <- OmegaHat[2, 1] <- rhohat * sqrt(omega1Hat * omega2Hat)

  c1 <- fit_step2$tau1[1]
  c2 <- fit_step2$tau2[1]

  zs1 <- z_scores1$V2[loci$idx_left[i_loci]:loci$idx_right[i_loci]][idx_include]
  zs2 <- z_scores2$V2[loci$idx_left[i_loci]:loci$idx_right[i_loci]][idx_include]

  fit_bg_ss1 <- XMAP::XMAP(simplify2array(list(R1)), cbind(zs1),
                           c(n1), K = K,
                           # Sigma1 = 0.1, Sigma2 = 0.1, rho = 0,
                           # Sig_E1 = 1, Sig_E2 = 1,
                           Omega = matrix(OmegaHat[1]),
                           Sig_E = c(c1),
                           # prior_weights = rep(1 / p, p),
                           tol = 1e-6,
                           maxIter = 200, estimate_residual_variance = F, estimate_prior_variance = T,
                           estimate_background_variance = F)

  ##################### XMAP bg with C=I #####################
  fit_CI <- XMAP(R1, R2, zs1, zs2, n1, n2, K = K,
                 Sigma1 = 0.1, Sigma2 = 0.1, rho = 0,
                 Sig_E1 = 1, Sig_E2 = 1,
                 Omega = OmegaHat,
                 prior_weights = rep(1 / p, p), tol = 1e-6,
                 maxIter = 200, estimate_residual_variance = F, estimate_prior_variance = T,
                 estimate_background_variance = F, initialize = F)
  cs1 <- get_CS(fit_CI, Xcorr = R1, coverage = 0.9, min_abs_corr = 0.1)
  cs2 <- get_CS(fit_CI, Xcorr = R2, coverage = 0.9, min_abs_corr = 0.1)
  cs <- cs1$cs[intersect(names(cs1$cs), names(cs2$cs))]
  pip <- get_pip(fit_CI$gamma)
  # plot_CS(pip, cs, Beta1[, i])
  cs_coverage[i, 1] <- mean(sapply(cs, function(x) any(idx_causal %in% x)))
  cs_power[i, 1] <- mean(idx_causal %in% unlist(cs))

  precision[i, 1] <- mean(which(pip > 0.9) %in% idx_causal)
  power[i, 1] <- mean(idx_causal %in% which(pip > 0.9))

  ##################### XMAP bg #####################
  fit_C <- XMAP(R1, R2, zs1, zs2, n1, n2, K = K,
                Sigma1 = 0.1, Sigma2 = 0.1, rho = 0,
                Sig_E1 = c1, Sig_E2 = c2,
                Omega = OmegaHat,
                prior_weights = rep(1 / p, p), tol = 1e-6,
                maxIter = 200, estimate_residual_variance = F, estimate_prior_variance = T,
                estimate_background_variance = F, initialize = F)
  fit_C <- XMAP::XMAP(simplify2array(list(R1,R2)), cbind(zs1,zs2),
                           c(n1,n2), K = K,
                           # Sigma1 = 0.1, Sigma2 = 0.1, rho = 0,
                           # Sig_E1 = 1, Sig_E2 = 1,
                           Omega = OmegaHat,
                           Sig_E = c(c1,c2),
                           # prior_weights = rep(1 / p, p),
                           tol = 1e-6,
                           maxIter = 200, estimate_residual_variance = F, estimate_prior_variance = T,
                           estimate_background_variance = F)
  cs1 <- get_CS(fit_C, Xcorr = R1, coverage = 0.9, min_abs_corr = 0.1)
  cs2 <- get_CS(fit_C, Xcorr = R2, coverage = 0.9, min_abs_corr = 0.1)
  cs <- cs1$cs[intersect(names(cs1$cs), names(cs2$cs))]
  pip <- get_pip(fit_C$gamma)
  # plot_CS(pip, cs, Beta1[, 1])
  cs_coverage[i, 2] <- mean(sapply(cs, function(x) any(idx_causal %in% x)))
  cs_power[i, 2] <- mean(idx_causal %in% unlist(cs))

  precision[i, 2] <- mean(which(pip > 0.9) %in% idx_causal)
  power[i, 2] <- mean(idx_causal %in% which(pip > 0.9))

  out_xmap_CI[[i]] <- fit_CI
  out_xmap_C[[i]] <- fit_C

  cat(i, "-th rep finished.\n")
}
saveRDS(out_xmap_CI, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/result/simu_inflate_pc1_", pc_var1, "_", pc_var2, "_loci_", i_loci, "_small_p", p, "_chr22_XMAP_CI_in_sample_n", n1, "_", n2,
                                   "_Ktrue", K_true, "_K", K, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                   "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".RDS"))

saveRDS(out_xmap_C, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/result/simu_inflate_pc1_", pc_var1, "_", pc_var2, "_loci_", i_loci, "_small_p", p, "_chr22_XMAP_C_in_sample_n", n1, "_", n2,
                                  "_Ktrue", K_true, "_K", K, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                  "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".RDS"))

saveRDS(out_ldscore, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/result/simu_inflate_pc1_", pc_var1, "_", pc_var2, "_loci_", i_loci, "_small_p", p, "_chr22_ldscore_in_sample_n", n1, "_", n2,
                                   "_Ktrue", K_true, "_K", K, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                   "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".RDS"))

