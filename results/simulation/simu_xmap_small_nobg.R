########################################################### pop1 ###########################################################
library(XPASS)
library(mvtnorm)
library(susieR)
library(data.table)
library(RhpcBLASctl)
set.seed(1)
blas_set_num_threads(15)
source("/import/home/share/mingxuan/fine_mapping/utils.R")
source("/import/home/share/mingxuan/fine_mapping/XMAP.R")
source("/home/share/mingxuan/taam/code/sldxr.R")
bim <- fread("/home/share/mingxuan/wegene_qc2/height_ukb_qc2_chr22.bim")
bp_chr <- bim$V4

nrep <- 50

# n1_all <- c(5000, 10000, 15000, 20000)
# n2_all <- c(10000, 15000, 20000)
n1 <- 20000
n2 <- 20000
K_true <- 3
K <- 5

h1_omega <- 0.005
h2_omega <- 0.005
rho_omega <- 0.8

enrich <- 50
enrich2 <- enrich*10
rho_sigma <- 0

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
R1 <- R1[idx_include,idx_include]
Beta1 <- data.matrix(fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_loci_",i_loci,"_small_p",p,"_Beta1_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                  "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".txt")))
R2 <- data.matrix(fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_R2_ukb_in_sample_ref_n", n2, "_chr22_loci", i_loci, ".txt")))
R2 <- R2[idx_include,idx_include]
Beta2 <- data.matrix(fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_loci_",i_loci,"_small_p",p,"_Beta2_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                  "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".txt")))


# output
cs_power <- cs_coverage <- power <- precision <- matrix(NA, nrep, 3)
out_xmap0 <- out_xmap_bg0 <- out_xmap_bg <- list()
for (i in 1:nrep) {
  idx_causal <- which(Beta1[, i] != 0)
  z_scores1 <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_zs1_nobg_loci_",i_loci,"_small_p",p,"_chr22_n", n1, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                            "_enrich", enrich2, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"))
  z_scores2 <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_zs2_nobg_loci_",i_loci,"_small_p",p,"_chr22_n", n2, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                            "_enrich", enrich2, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"))
  susie_plot(z_scores1$V2,"z",b=Beta1[, i])
  susie_plot(z_scores2$V2,"z",b=Beta2[, i])

  # ldsc estimate parameters
  idx1 <- which(z_scores1$V2^2 < 30 & z_scores2$V2^2 < 30)
  ld1 <- colSums(R1^2)
  ld2 <- colSums(R2^2)
  ld1 <- ld1 - (1 - ld1) / (n1 - 2)
  ld2 <- ld2 - (1 - ld2) / (n2 - 2)
  ld12 <- colSums(R1 * R2)

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

  ##################### XMAP bg with rho=0 #####################
  fit_bg_ss <- XMAP(R1, R2, z_scores1$V2, z_scores2$V2, n1, n2, K = K,
                    Sigma1 = 0.1, Sigma2 = 0.1, rho = 0,
                    Sig_E1 = 0.9, Sig_E2 = 0.9,
                    Omega = OmegaHat0,
                    prior_weights = rep(1 / p, p), tol = 1e-6,
                    maxIter = 200, estimate_residual_variance = T, estimate_prior_variance = T,
                    estimate_background_variance = F, initialize = F)
  cs1 <- get_CS(fit_bg_ss, Xcorr = R1, coverage = 0.9, min_abs_corr = 0.1)
  cs2 <- get_CS(fit_bg_ss, Xcorr = R2, coverage = 0.9, min_abs_corr = 0.1)
  cs <- cs1$cs[intersect(names(cs1$cs), names(cs2$cs))]
  pip <- get_pip(fit_bg_ss$gamma)
  # plot_CS(pip, cs, Beta1[, i])
  cs_coverage[i, 1] <- mean(sapply(cs, function(x) any(idx_causal %in% x)))
  cs_power[i, 1] <- mean(idx_causal %in% unlist(cs))

  precision[i, 1] <- mean(which(pip > 0.9) %in% idx_causal)
  power[i, 1] <- mean(idx_causal %in% which(pip > 0.9))

  ##################### XMAP bg #####################
  fit_bg_ss1 <- XMAP(R1, R2, z_scores1$V2, z_scores2$V2, n1, n2, K = K,
                     Sigma1 = 0.1, Sigma2 = 0.1, rho = 0,
                     Sig_E1 = 1, Sig_E2 = 1,
                     Omega = OmegaHat,
                     prior_weights = rep(1 / p, p), tol = 1e-6,
                     maxIter = 200, estimate_residual_variance = T, estimate_prior_variance = T,
                     estimate_background_variance = F, initialize = F)
  cs1 <- get_CS(fit_bg_ss1, Xcorr = R1, coverage = 0.9, min_abs_corr = 0.1)
  cs2 <- get_CS(fit_bg_ss1, Xcorr = R2, coverage = 0.9, min_abs_corr = 0.1)
  cs <- cs1$cs[intersect(names(cs1$cs), names(cs2$cs))]
  pip <- get_pip(fit_bg_ss1$gamma)
  # plot_CS(pip, cs, Beta1[, 1])
  cs_coverage[i, 2] <- mean(sapply(cs, function(x) any(idx_causal %in% x)))
  cs_power[i, 2] <- mean(idx_causal %in% unlist(cs))

  precision[i, 2] <- mean(which(pip > 0.9) %in% idx_causal)
  power[i, 2] <- mean(idx_causal %in% which(pip > 0.9))

  ##################### XMAP no bg #####################
  fit_ss <- XMAP(R1, R2, z_scores1$V2, z_scores2$V2, n1, n2, K = K,
                 Sigma1 = 0.1, Sigma2 = 0.1, rho = 0,
                 Sig_E1 = 1, Sig_E2 = 1,
                 include_background = F,
                 prior_weights = rep(1 / p, p), tol = 1e-6,
                 maxIter = 200, estimate_residual_variance = F, estimate_prior_variance = T,
                 estimate_background_variance = F, initialize = F)
  cs1 <- get_CS(fit_ss, Xcorr = R1, coverage = 0.9, min_abs_corr = 0.1)
  cs2 <- get_CS(fit_ss, Xcorr = R2, coverage = 0.9, min_abs_corr = 0.1)
  cs <- cs1$cs[intersect(names(cs1$cs), names(cs2$cs))]
  pip <- get_pip(fit_ss$gamma)
  cs_coverage[i, 3] <- mean(sapply(cs, function(x) any(idx_causal %in% x)))
  cs_power[i, 3] <- mean(idx_causal %in% unlist(cs))

  precision[i, 3] <- mean(which(pip > 0.9) %in% idx_causal)
  power[i, 3] <- mean(idx_causal %in% which(pip > 0.9))

  out_xmap0[[i]] <- fit_ss
  out_xmap_bg[[i]] <- fit_bg_ss1
  out_xmap_bg0[[i]] <- fit_bg_ss


  cat(i, "-th rep finished.\n")
}
saveRDS(out_xmap0, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/result/simu_nobg_loci_",i_loci,"_small_p",p,"_chr22_XMAP0_in_sample_n", n1, "_", n2,
                                        "_Ktrue", K_true, "_K", K, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                        "_enrich", enrich2, "_rhoSigma", rho_sigma, "_nrep", nrep, ".RDS"))

saveRDS(out_xmap_bg, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/result/simu_nobg_loci_",i_loci,"_small_p",p,"_chr22_XMAP_bg_in_sample_n", n1, "_", n2,
                                          "_Ktrue", K_true, "_K", K, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                          "_enrich", enrich2, "_rhoSigma", rho_sigma, "_nrep", nrep, ".RDS"))

saveRDS(out_xmap_bg0, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/result/simu_nobg_loci_",i_loci,"_small_p",p,"_chr22_XMAP_bg0_in_sample_n", n1, "_", n2,
                                          "_Ktrue", K_true, "_K", K, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                          "_enrich", enrich2, "_rhoSigma", rho_sigma, "_nrep", nrep, ".RDS"))


write.table(cs_power, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/result/simu_nobg_loci_",i_loci,"_small_p",p,"_chr22_XMAP_in_sample_cs_power_n", n1, "_", "_", n2,
                                    "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                    "_enrich", enrich2, "_rhoSigma", rho_sigma, "_nrep", nrep, "_cvg0.9_minCor0.1.txt"),
            col.names = F, row.names = F, quote = F, sep = "\t")
write.table(cs_coverage, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/result/simu_nobg_loci_",i_loci,"_small_p",p,"_chr22_XMAP_in_sample_cs_coverage_n", n1, "_", n2,
                                       "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                       "_enrich", enrich2, "_rhoSigma", rho_sigma, "_nrep", nrep, "_cvg0.9_minCor0.1.txt"),
            col.names = F, row.names = F, quote = F, sep = "\t")


write.table(power, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/result/simu_nobg_loci_",i_loci,"_small_p",p,"_chr22_XMAP_in_sample_power_n", n1, "_", n2,
                                 "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                 "_enrich", enrich2, "_rhoSigma", rho_sigma, "_nrep", nrep, "_pip0.9.txt"),
            col.names = F, row.names = F, quote = F, sep = "\t")
write.table(precision, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/result/simu_nobg_loci_",i_loci,"_small_p",p,"_chr22_XMAP_in_sample_precision_n", n1, "_", n2,
                                     "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                     "_enrich", enrich2, "_rhoSigma", rho_sigma, "_nrep", nrep, "_pip0.9.txt"),
            col.names = F, row.names = F, quote = F, sep = "\t")

