########################################################### pop1 ###########################################################
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
source("/import/home/share/mingxuan/taam/code/sldxr.R")
bim <- fread("/import/home/share/mingxuan/wegene_qc2/height_ukb_qc2_chr22.bim")
bp_chr <- bim$V4

nrep <- 50

# n1_all <- c(5000, 10000, 15000, 20000)
# n2_all <- c(10000, 15000, 20000)
n1 <- 5000
n2 <- 20000
K_true <- 3
K <- 10

h1_omega <- 0.005
h2_omega <- 0.005
rho_omega <- 0.8

enrich <- 50
rho_sigma <- 0

loci <- fread(paste0("/import/home/share/mingxuan/fine_mapping/loci_chr22_1mb.txt"))
i_loci <- 29
idx_include <- (501:1000) + loci$idx_start[i_loci] - 1
p <- length(idx_include)


# per-snp variance for bg
omega1 <- h1_omega / p
omega2 <- h2_omega / p

Omega <- matrix(c(omega1, rho_omega * sqrt(omega1 * omega2), rho_omega * sqrt(omega1 * omega2), omega2), 2, 2)

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

ldscore <- fread(paste0("/import/home/share/mingxuan/fine_mapping/LD_mat/LDscore_UKB_wegene_chr22.txt"))

# output
cs_power <- cs_coverage <- power <- precision <- matrix(NA, nrep, 3)
out_susie1 <- out_susie2 <- out_xmap0 <- out_xmap <- out_xmap1 <- out_xmap2 <- out_xmap1_local <- out_xmap2_local <- out_xmap_local <- out_xmap_true <- list()
Omega_all <- Omega_raw_all <- C_all <- list()
for (i in 1:nrep) {
  idx_causal <- which(Beta1[, i] != 0)
  z_scores1 <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_zs1_wholeCHR_loci_", i_loci, "_small_p", p, "_chr22_n", n1, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                            "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"))
  z_scores2 <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_zs2_wholeCHR_loci_", i_loci, "_small_p", p, "_chr22_n", n2, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                            "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"))


  # ldsc estimate parameters
  idx1 <- 1:nrow(z_scores1) # which(z_scores1$V2^2 < 30 & z_scores2$V2^2 < 30)
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

  Omega_raw_all[[i]] <- matrix(c(fit_step2$tau1[2], fit_step2$theta[2], fit_step2$theta[2], fit_step2$tau2[2]), 2, 2)
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

  Omega_all[[i]] <- OmegaHat
  C_all[[i]] <- c(fit_step2$tau1[1], fit_step2$tau2[1])

  zs1 <- z_scores1$V2[loci$idx_left[i_loci]:loci$idx_right[i_loci]][idx_include]
  zs2 <- z_scores2$V2[loci$idx_left[i_loci]:loci$idx_right[i_loci]][idx_include]




  ##################### XMAP pop1 fix bg #####################
  fit_xmap1 <- XMAP::XMAP(simplify2array(list(R1)), cbind(zs1),
                          c(n1), K = K,
                          # Sigma1 = 0.1, Sigma2 = 0.1, rho = 0,
                          # Sig_E1 = 1, Sig_E2 = 1,
                          Omega = matrix(OmegaHat[1, 1]),
                          Sig_E = c(1),
                          # prior_weights = rep(1 / p, p),
                          tol = 1e-8,
                          maxIter = 300, estimate_residual_variance = F, estimate_prior_variance = T,
                          estimate_background_variance = F)
  ##################### XMAP pop1 update bg #####################
  fit_xmap1_local <- XMAP::XMAP(simplify2array(list(R1)), cbind(zs1),
                                c(n1), K = K,
                                # Sigma1 = 0.1, Sigma2 = 0.1, rho = 0,
                                # Sig_E1 = 1, Sig_E2 = 1,
                                Omega = matrix(OmegaHat[1, 1]),
                                Sig_E = c(1),
                                # prior_weights = rep(1 / p, p),
                                tol = 1e-8,
                                maxIter = 300, estimate_residual_variance = F, estimate_prior_variance = T,
                                estimate_background_variance = T)
  ##################### XMAP pop2 fix bg #####################
  fit_xmap2 <- XMAP::XMAP(simplify2array(list(R2)), cbind(zs2),
                          c(n2), K = K,
                          # Sigma1 = 0.1, Sigma2 = 0.1, rho = 0,
                          # Sig_E1 = 1, Sig_E2 = 1,
                          Omega = matrix(OmegaHat[2, 2]),
                          Sig_E = c(1),
                          # prior_weights = rep(1 / p, p),
                          tol = 1e-8,
                          maxIter = 300, estimate_residual_variance = F, estimate_prior_variance = T,
                          estimate_background_variance = F)
  ##################### XMAP pop2 update bg #####################
  fit_xmap2_local <- XMAP::XMAP(simplify2array(list(R2)), cbind(zs2),
                                c(n2), K = K,
                                # Sigma1 = 0.1, Sigma2 = 0.1, rho = 0,
                                # Sig_E1 = 1, Sig_E2 = 1,
                                Omega = matrix(OmegaHat[2, 2]),
                                Sig_E = c(1),
                                # prior_weights = rep(1 / p, p),
                                tol = 1e-8,
                                maxIter = 300, estimate_residual_variance = F, estimate_prior_variance = T,
                                estimate_background_variance = T)
  ##################### XMAP fix bg #####################
  fit_xmap <- XMAP::XMAP(simplify2array(list(R1, R2)), cbind(zs1, zs2),
                         c(n1, n2), K = K,
                         # Sigma1 = 0.1, Sigma2 = 0.1, rho = 0,
                         # Sig_E1 = 1, Sig_E2 = 1,
                         Omega = OmegaHat,
                         Sig_E = c(1, 1),
                         # prior_weights = rep(1 / p, p),
                         tol = 1e-8,
                         maxIter = 300, estimate_residual_variance = F, estimate_prior_variance = T,
                         estimate_background_variance = F)
  ##################### XMAP update bg #####################
  fit_xmap_local <- XMAP::XMAP(simplify2array(list(R1, R2)), cbind(zs1, zs2),
                               c(n1, n2), K = K,
                               # Sigma1 = 0.1, Sigma2 = 0.1, rho = 0,
                               # Sig_E1 = 1, Sig_E2 = 1,
                               Omega = OmegaHat,
                               Sig_E = c(1, 1),
                               # prior_weights = rep(1 / p, p),
                               tol = 1e-8,
                               maxIter = 300, estimate_residual_variance = F, estimate_prior_variance = T,
                               estimate_background_variance = T)

  ##################### XMAP true bg #####################
  fit_xmap_true <- XMAP::XMAP(simplify2array(list(R1, R2)), cbind(zs1, zs2),
                              c(n1, n2), K = K,
                              # Sigma1 = 0.1, Sigma2 = 0.1, rho = 0,
                              # Sig_E1 = 1, Sig_E2 = 1,
                              Omega = Omega,
                              Sig_E = c(1, 1),
                              # prior_weights = rep(1 / p, p),
                              tol = 1e-8,
                              maxIter = 300, estimate_residual_variance = F, estimate_prior_variance = T,
                              estimate_background_variance = T)

  ##################### XMAP no bg #####################
  fit_xmap0 <- XMAP::XMAP(simplify2array(list(R1, R2)), cbind(zs1, zs2),
                          c(n1, n2), K = K,
                          Omega = matrix(c(1e-30, 0, 0, 1e-30), 2, 2),
                          Sig_E = c(1, 1),
                          # prior_weights = rep(1 / p, p),
                          tol = 1e-8,
                          maxIter = 300, estimate_residual_variance = F, estimate_prior_variance = T,
                          estimate_background_variance = F)


  out_xmap_true[[i]] <- fit_xmap_true
  out_xmap0[[i]] <- fit_xmap0
  out_xmap1[[i]] <- fit_xmap1
  out_xmap1_local[[i]] <- fit_xmap1_local
  out_xmap2[[i]] <- fit_xmap2
  out_xmap2_local[[i]] <- fit_xmap2_local
  out_xmap[[i]] <- fit_xmap
  out_xmap_local[[i]] <- fit_xmap_local


  cat(i, "-th rep finished.\n")
}
saveRDS(Omega_all, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/result/simu_WholeCHR_loci_", i_loci, "_small_p", p, "_chr22_OmegaHat_in_sample_n", n1,
                                 "_Ktrue", K_true, "_K", K, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                 "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".RDS"))

saveRDS(C_all, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/result/simu_WholeCHR_loci_", i_loci, "_small_p", p, "_chr22_CHat_in_sample_n", n1,
                             "_Ktrue", K_true, "_K", K, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                             "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".RDS"))



saveRDS(out_xmap1, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/result/simu_WholeCHR_loci_", i_loci, "_small_p", p, "_chr22_XMAP_bg_single1_in_sample_n", n1, "_", n2,
                                 "_Ktrue", K_true, "_K", K, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                 "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".RDS"))
saveRDS(out_xmap2, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/result/simu_WholeCHR_loci_", i_loci, "_small_p", p, "_chr22_XMAP_bg_single2_in_sample_n", n1, "_", n2,
                                 "_Ktrue", K_true, "_K", K, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                 "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".RDS"))

saveRDS(out_xmap1_local, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/result/simu_WholeCHR_loci_", i_loci, "_small_p", p, "_chr22_XMAP_bg_single1_local_in_sample_n", n1, "_", n2,
                                       "_Ktrue", K_true, "_K", K, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                       "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".RDS"))
saveRDS(out_xmap2_local, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/result/simu_WholeCHR_loci_", i_loci, "_small_p", p, "_chr22_XMAP_bg_single2_local_in_sample_n", n1, "_", n2,
                                       "_Ktrue", K_true, "_K", K, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                       "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".RDS"))

saveRDS(out_xmap, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/result/simu_WholeCHR_loci_", i_loci, "_small_p", p, "_chr22_XMAP_bg_in_sample_n", n1, "_", n2,
                                "_Ktrue", K_true, "_K", K, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".RDS"))

saveRDS(out_xmap_local, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/result/simu_WholeCHR_loci_", i_loci, "_small_p", p, "_chr22_XMAP_bg_local_in_sample_n", n1, "_", n2,
                                      "_Ktrue", K_true, "_K", K, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                      "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".RDS"))

saveRDS(out_xmap0, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/result/simu_WholeCHR_loci_", i_loci, "_small_p", p, "_chr22_XMAP_bg0_in_sample_n", n1, "_", n2,
                                 "_Ktrue", K_true, "_K", K, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                 "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".RDS"))

saveRDS(out_xmap_true, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/result/simu_WholeCHR_loci_", i_loci, "_small_p", p, "_chr22_XMAP_bg_true_in_sample_n", n1, "_", n2,
                                     "_Ktrue", K_true, "_K", K, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                     "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".RDS"))
