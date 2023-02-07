#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

library(data.table)
library(susieR)
library(XMAP)


sumstat_UKB <- fread("/home/share/sumstats/format/Height_allSNPs_UKBNealLab_summary_format.txt")
sumstat_WG <- fread("/home/share/sumstats/format/height_Chinese_allSNPs_summary_format.txt")
sumstat_sibship <- fread("/home/share/sumstats/format/height_EUR_sibship_summary_format.txt")
# sumstat_GIANT <- fread("/home/share/sumstats/format/Height_allSNPs_UKBNealLab_summary_format.txt")
sumstat_BBJ <- fread("/home/share/sumstats/format/height_BBJ2020_beta_se_summary_MTAG_format.txt")

ldscore <- data.frame()
for (chr in 1:22) {
  ldscore_chr <- fread(paste0("/home/share/mingxuan/UKB_geno_finemap/LDscore_eas_brit_afr_chr", chr, ".txt"))
  ldscore <- rbind(ldscore, ldscore_chr)
  cat("CHR", chr, "\n")
}
# remove ambiguous SNPs
idx_amb <- which(ldscore$allele1 == comple(ldscore$allele2))
ldscore <- ldscore[-idx_amb,]


# overlap SNPs
snps <- Reduce(intersect, list(ldscore$rsid, sumstat_UKB$SNP, sumstat_WG$SNP, sumstat_sibship$SNP, sumstat_BBJ$SNP))
sumstat_UKB <- sumstat_UKB[match(snps, sumstat_UKB$SNP),]
sumstat_WG <- sumstat_WG[match(snps, sumstat_WG$SNP),]
sumstat_sibship <- sumstat_sibship[match(snps, sumstat_sibship$SNP),]
sumstat_BBJ <- sumstat_BBJ[match(snps, sumstat_BBJ$SNP),]
ldscore <- ldscore[match(snps, ldscore$rsid),]

# write.table(snps,file="/home/share/mingxuan/fine_mapping/analysis/snps_height_UKB_sibship_WG_BBJ.txt",col.names = F,row.names = F,quote=F)

# flip alleles
z_ukb <- sumstat_UKB$Z
# z_eur <- sumstat_EUR$beta / sumstat_EUR$se
z_wg <- sumstat_WG$Z
z_sibship <- sumstat_sibship$Z
z_bbj <- sumstat_BBJ$Z

idx_flip <- which(sumstat_UKB$A1 != ldscore$allele1 & sumstat_UKB$A1 != comple(ldscore$allele1))
z_ukb[idx_flip] <- -z_ukb[idx_flip]

idx_flip <- which(sumstat_WG$A1 != ldscore$allele1 & sumstat_WG$A1 != comple(ldscore$allele1))
z_wg[idx_flip] <- -z_wg[idx_flip]

idx_flip <- which(sumstat_sibship$A1 != ldscore$allele1 & sumstat_sibship$A1 != comple(ldscore$allele1))
z_sibship[idx_flip] <- -z_sibship[idx_flip]

idx_flip <- which(sumstat_BBJ$A1 != ldscore$allele1 & sumstat_BBJ$A1 != comple(ldscore$allele1))
z_bbj[idx_flip] <- -z_bbj[idx_flip]


######################### NealLab UKB GWAS + Chinese GWAS #########################
idx1 <- which(z_ukb^2 < 30 & z_wg^2 < 30)

ld_eas_w <- 1 / sapply(ldscore$EAS, function(x) max(x, 1))
ld_eur_w <- 1 / sapply(ldscore$EUR, function(x) max(x, 1))
fit_step1 <- estimate_gc(data.frame(Z = z_wg[idx1], N = sumstat_WG$N[idx1]), data.frame(Z = z_ukb[idx1], N = sumstat_UKB$N[idx1]),
                         ldscore$EAS[idx1], ldscore$EUR[idx1], ldscore$EUR_EAS[idx1],
                         reg_w1 = ld_eas_w[idx1], reg_w2 = ld_eur_w[idx1], reg_wx = sqrt(ld_eas_w[idx1] * ld_eur_w[idx1]),
                         constrain_intercept = F)


fit_step2 <- estimate_gc(data.frame(Z = z_wg, N = sumstat_WG$N), data.frame(Z = z_ukb, N = sumstat_UKB$N),
                         ldscore$EAS, ldscore$EUR, ldscore$EUR_EAS,
                         reg_w1 = ld_eas_w, reg_w2 = ld_eur_w, reg_wx = sqrt(ld_eas_w * ld_eur_w),
                         constrain_intercept = T, fit_step1$tau1$coefs[1], fit_step1$tau2$coefs[1], fit_step1$theta$coefs[1])

######################### LDSC output #########################
# [1] 1.069005e+00 8.448969e-08
# $tau2
# [1] 1.658928e+00 7.596966e-08
# $theta
# [1] 9.838635e-02 6.990783e-08
###################################################################


######################### NealLab UKB GWAS + BBJ GWAS #########################
idx1 <- which(z_ukb^2 < 30 & z_bbj^2 < 30)

ld_eas_w <- 1 / sapply(ldscore$EAS, function(x) max(x, 1))
ld_eur_w <- 1 / sapply(ldscore$EUR, function(x) max(x, 1))
fit_step1 <- estimate_gc(data.frame(Z = z_bbj[idx1], N = sumstat_BBJ$N[idx1]), data.frame(Z = z_ukb[idx1], N = sumstat_UKB$N[idx1]),
                         ldscore$EAS[idx1], ldscore$EUR[idx1], ldscore$EUR_EAS[idx1],
                         reg_w1 = ld_eas_w[idx1], reg_w2 = ld_eur_w[idx1], reg_wx = sqrt(ld_eas_w[idx1] * ld_eur_w[idx1]),
                         constrain_intercept = F)


fit_step2 <- estimate_gc(data.frame(Z = z_bbj, N = sumstat_BBJ$N), data.frame(Z = z_ukb, N = sumstat_UKB$N),
                         ldscore$EAS, ldscore$EUR, ldscore$EUR_EAS,
                         reg_w1 = ld_eas_w, reg_w2 = ld_eur_w, reg_wx = sqrt(ld_eas_w * ld_eur_w),
                         constrain_intercept = T, fit_step1$tau1$coefs[1], fit_step1$tau2$coefs[1], fit_step1$theta$coefs[1])

######################### LDSC output #########################
# $tau1
# [1] 1.376723e+00 8.013839e-08
# $tau2
# [1] 1.656730e+00 7.607146e-08
# $theta
# [1] 1.806504e-01 6.784428e-08
###################################################################


######################### Sibship GWAS + WG GWAS #########################
idx1 <- which(z_sibship^2 < 30 & z_wg^2 < 30)

ld_eas_w <- 1 / sapply(ldscore$EAS, function(x) max(x, 1))
ld_eur_w <- 1 / sapply(ldscore$EUR, function(x) max(x, 1))
fit_step1 <- estimate_gc(data.frame(Z = z_wg[idx1], N = sumstat_WG$N[idx1]), data.frame(Z = z_sibship[idx1], N = sumstat_sibship$N[idx1]),
                         ldscore$EAS[idx1], ldscore$EUR[idx1], ldscore$EUR_EAS[idx1],
                         reg_w1 = ld_eas_w[idx1], reg_w2 = ld_eur_w[idx1], reg_wx = sqrt(ld_eas_w[idx1] * ld_eur_w[idx1]),
                         constrain_intercept = F)


fit_step2 <- estimate_gc(data.frame(Z = z_wg, N = sumstat_WG$N), data.frame(Z = z_sibship, N = sumstat_sibship$N),
                         ldscore$EAS, ldscore$EUR, ldscore$EUR_EAS,
                         reg_w1 = ld_eas_w, reg_w2 = ld_eur_w, reg_wx = sqrt(ld_eas_w * ld_eur_w),
                         constrain_intercept = T, fit_step1$tau1$coefs[1], fit_step1$tau2$coefs[1], fit_step1$theta$coefs[1])

######################### LDSC output #########################
# $tau1
# [1] 1.073431e+00 8.267105e-08
# $tau2
# [1] 1.114412e+00 5.436051e-08
# $theta
# [1] 4.874585e-02 5.674659e-08
###################################################################


######################### Sibship GWAS + BBJ GWAS #########################
idx1 <- which(z_sibship^2 < 30 & z_bbj^2 < 30)

ld_eas_w <- 1 / sapply(ldscore$EAS, function(x) max(x, 1))
ld_eur_w <- 1 / sapply(ldscore$EUR, function(x) max(x, 1))
fit_step1 <- estimate_gc(data.frame(Z = z_bbj[idx1], N = sumstat_BBJ$N[idx1]), data.frame(Z = z_sibship[idx1], N = sumstat_sibship$N[idx1]),
                         ldscore$EAS[idx1], ldscore$EUR[idx1], ldscore$EUR_EAS[idx1],
                         reg_w1 = ld_eas_w[idx1], reg_w2 = ld_eur_w[idx1], reg_wx = sqrt(ld_eas_w[idx1] * ld_eur_w[idx1]),
                         constrain_intercept = F)


fit_step2 <- estimate_gc(data.frame(Z = z_bbj, N = sumstat_BBJ$N), data.frame(Z = z_sibship, N = sumstat_sibship$N),
                         ldscore$EAS, ldscore$EUR, ldscore$EUR_EAS,
                         reg_w1 = ld_eas_w, reg_w2 = ld_eur_w, reg_wx = sqrt(ld_eas_w * ld_eur_w),
                         constrain_intercept = T, fit_step1$tau1$coefs[1], fit_step1$tau2$coefs[1], fit_step1$theta$coefs[1])

######################### LDSC output #########################
# $tau1
# [1] 1.394916e+00 7.809809e-08
# $tau2
# [1] 1.116069e+00 5.411153e-08
# $theta
# [1] 9.205534e-02 5.591833e-08
###################################################################


library(data.table)
library(RhpcBLASctl)
library(susieR)
library(Matrix)
blas_set_num_threads(40)
# library(susieR)
source("/import/home/share/mingxuan/fine_mapping/utils.R")
source("/import/home/share/mingxuan/fine_mapping/XMAP.R")
source("/home/share/mingxuan/taam/code/sldxr.R")
sumstat_UKB <- fread("/import/home/share/sumstats/format/Height_allSNPs_UKBNealLab_summary_format.txt")
sumstat_WG <- fread("/import/home/share/sumstats/format/height_Chinese_allSNPs_summary_format.txt")
sumstat_sibship <- fread("/import/home/share/sumstats/format/height_EUR_sibship_summary_format.txt")
# sumstat_GIANT <- fread("/import/home/share/sumstats/format/Height_allSNPs_UKBNealLab_summary_format.txt")
sumstat_BBJ <- fread("/import/home/share/sumstats/format/height_BBJ2020_beta_se_summary_MTAG_format.txt")

snps <- fread("/import/home/share/mingxuan/fine_mapping/analysis/snps_height_UKB_sibship_WG_BBJ.txt", header = F)$V1
sumstat_UKB <- sumstat_UKB[match(snps, sumstat_UKB$SNP),]
sumstat_WG <- sumstat_WG[match(snps, sumstat_WG$SNP),]
sumstat_sibship <- sumstat_sibship[match(snps, sumstat_sibship$SNP),]
sumstat_BBJ <- sumstat_BBJ[match(snps, sumstat_BBJ$SNP),]

K <- 15

# eur_data <- "UKB"
eur_data <- "Sibship"

# eas_data <- "WG"
eas_data <- "BBJ"

if (eur_data == "UKB" & eas_data == "WG") {
  ################################# NealLab UKB GWAS + Chinese GWAS ##################################
  sumstat_EUR <- sumstat_UKB
  sumstat_EAS <- sumstat_WG
  # $tau1
  # [1] 1.069005e+00 8.448969e-08
  # $tau2
  # [1] 1.658928e+00 7.596966e-08
  # $theta
  # [1] 9.838635e-02 6.990783e-08
  c1 <- 1.069005
  c2 <- 1.658928
  Omega1 <- 8.448969e-08
  Omega2 <- 7.596966e-08
  Omega12 <- 6.990783e-08
  OmegaHat <- OmegaHat0 <- diag(c(Omega1, Omega2))
  OmegaHat[1, 2] <- OmegaHat[2, 1] <- Omega12
  #####################################################################################################
}

if (eur_data == "UKB" & eas_data == "BBJ") {
  ################################# NealLab UKB GWAS + BBJ GWAS ##################################
  sumstat_EUR <- sumstat_UKB
  sumstat_EAS <- sumstat_BBJ
  # $tau1
  # [1] 1.376723e+00 8.013839e-08
  # $tau2
  # [1] 1.656730e+00 7.607146e-08
  # $theta
  # [1] 1.806504e-01 6.784428e-08
  c1 <- 1.376723
  c2 <- 1.656730
  Omega1 <- 8.013839e-08
  Omega2 <- 7.607146e-08
  Omega12 <- 6.784428e-08
  OmegaHat <- OmegaHat0 <- diag(c(Omega1, Omega2))
  OmegaHat[1, 2] <- OmegaHat[2, 1] <- Omega12
  #####################################################################################################
}


if (eur_data == "Sibship" & eas_data == "WG") {
  ################################# Sibship GWAS + WG GWAS ##################################
  sumstat_EUR <- sumstat_sibship
  sumstat_EAS <- sumstat_WG
  # $tau1
  # [1] 1.073431e+00 8.267105e-08
  # $tau2
  # [1] 1.114412e+00 5.436051e-08
  # $theta
  # [1] 4.874585e-02 5.674659e-08
  c1 <- 1.073431
  c2 <- 1.114412
  Omega1 <- 8.267105e-08
  Omega2 <- 5.436051e-08
  Omega12 <- 5.674659e-08
  OmegaHat <- OmegaHat0 <- diag(c(Omega1, Omega2))
  OmegaHat[1, 2] <- OmegaHat[2, 1] <- Omega12
  #####################################################################################################
}


if (eur_data == "Sibship" & eas_data == "BBJ") {
  ################################# Sibship GWAS + BBJ GWAS ##################################
  sumstat_EUR <- sumstat_sibship
  sumstat_EAS <- sumstat_BBJ
  # $tau1
  # [1] 1.394916e+00 7.809809e-08
  # $tau2
  # [1] 1.116069e+00 5.411153e-08
  # $theta
  # [1] 9.205534e-02 5.591833e-08
  c1 <- 1.394916
  c2 <- 1.116069
  Omega1 <- 7.809809e-08
  Omega2 <- 5.411153e-08
  Omega12 <- 5.591833e-08
  OmegaHat <- OmegaHat0 <- diag(c(Omega1, Omega2))
  OmegaHat[1, 2] <- OmegaHat[2, 1] <- Omega12
  #####################################################################################################
}

chr_all <- 1:22 #as.numeric(args[1])
for (chr in chr_all) {
  loci <- fread(paste0("/import/home/share/mingxuan/UKB_geno_finemap/chr", chr, "_eas_brit_afr.loci"))
  for (i in 1:nrow(loci)) {

    info <- fread(paste0("/import/home/share/mingxuan/UKB_geno_finemap/chr", chr, "_", loci$left[i], "_", loci$right[i], "_eas.info"))
    info_eur_afr <- fread(paste0("/import/home/share/mingxuan/UKB_geno_finemap/chr", chr, "_", loci$left[i], "_", loci$right[i], ".info"))
    # remove ambiguous SNPs
    idx_amb <- which(info$allele1 == comple(info$allele2))
    info <- info[-idx_amb,]


    snps <- Reduce(intersect, list(info$rsid, sumstat_EAS$SNP, sumstat_EUR$SNP))

    if (length(snps) < 10) next
    # overlap SNPs
    sumstat_EAS_i <- sumstat_EAS[match(snps, sumstat_EAS$SNP),]
    sumstat_EUR_i <- sumstat_EUR[match(snps, sumstat_EUR$SNP),]
    # if (any(2 * pnorm(abs(sumstat_EAS_i$beta / sumstat_EAS_i$se), lower.tail = F) <= 5e-8) | all(2 * pnorm(abs(sumstat_EUR_i$Z), lower.tail = F) <= 5e-8)) {
    #   cat(i,"-th loci has no significant snps.\n")
    #   next
    # }
    R_eas <- readMM(paste0("/import/home/share/mingxuan/UKB_geno_finemap/chr", chr, "_", loci$left[i], "_", loci$right[i], "_eas.mtx.gz"))
    R_eas <- as.matrix(R_eas + t(R_eas))
    R_eas <- R_eas[-idx_amb, -idx_amb]
    idx_eas <- match(snps, info$rsid)
    R_eas <- R_eas[idx_eas, idx_eas]
    info <- info[idx_eas,]

    R_brit <- readMM(paste0("/import/home/share/mingxuan/UKB_geno_finemap/chr", chr, "_", loci$left[i], "_", loci$right[i], "_brit.mtx.gz"))
    R_brit <- as.matrix(R_brit + t(R_brit))
    idx_brit <- match(snps, info_eur_afr$rsid)
    R_brit <- R_brit[idx_brit, idx_brit]

    # remove sample size outliers
    q_EUR <- quantile(sumstat_EUR_i$N)
    iqr_EUR <- IQR(sumstat_EUR_i$N)
    idx_outlier_EUR <- which(sumstat_EUR_i$N < q_EUR[2] - iqr_EUR | sumstat_EUR_i$N > q_EUR[3] + iqr_EUR)

    q_EAS <- quantile(sumstat_EAS_i$N)
    iqr_EAS <- IQR(sumstat_EAS_i$N)
    idx_outlier_EAS <- which(sumstat_EAS_i$N < q_EAS[2] - iqr_EAS | sumstat_EAS_i$N > q_EAS[3] + iqr_EAS)

    idx_outlier <- unique(c(idx_outlier_EUR, idx_outlier_EAS))
    if (length(idx_outlier) > 0) {

      snps <- snps[-idx_outlier]
      sumstat_EAS_i <- sumstat_EAS_i[-idx_outlier,]
      sumstat_EUR_i <- sumstat_EUR_i[-idx_outlier,]
      info <- info[-idx_outlier,]
      R_eas <- R_eas[-idx_outlier, -idx_outlier]
      R_brit <- R_brit[-idx_outlier, -idx_outlier]
    }

    # flip alleles
    z_eas <- sumstat_EAS_i$Z
    z_eur <- sumstat_EUR_i$Z

    idx_flip <- which(sumstat_EAS_i$A1 != info$allele1 & sumstat_EAS_i$A1 != comple(info$allele1))
    z_eas[idx_flip] <- -z_eas[idx_flip]

    idx_flip <- which(sumstat_EUR_i$A1 != info$allele1 & sumstat_EUR_i$A1 != comple(info$allele1))
    z_eur[idx_flip] <- -z_eur[idx_flip]


    susie_eas <- susie_rss(z_eas, R_eas,
                           L = K,
                           verbose = TRUE,
                           coverage = 0.9,
                           min_abs_corr = 0.1,
                           check_prior = F)
    saveRDS(susie_eas, file = paste0("/import/home/share/mingxuan/fine_mapping/analysis/results/height_susie_K",K,"_out_", eas_data, "_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))
    
    susie_eur <- susie_rss(z_eur, R_brit,
                           L = K,
                           verbose = TRUE,
                           coverage = 0.9,
                           min_abs_corr = 0.1,
                           check_prior = F)
    saveRDS(susie_eur, file = paste0("/import/home/share/mingxuan/fine_mapping/analysis/results/height_susie_K",K,"_out_", eur_data, "_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))
    # # susie_plot(susie_eur, y = "PIP",  main = "SuSiE-UKB_EUR")
    # # susie_plot(z_eur, y = "z", main = "GWAS-UKB_EUR")
    #
    #
    # # full XMAP model
    xmap_C <- XMAP::XMAP(simplify2array(list(R_eas, R_brit)), cbind(z_eas, z_eur), c(median(sumstat_EAS_i$N), median(sumstat_EUR_i$N)),
                         K = K, Omega = OmegaHat, Sig_E = c(c1, c2), tol = 1e-6,
                         maxIter = 200, estimate_residual_variance = F, estimate_prior_variance = T,
                         estimate_background_variance = F)
    saveRDS(xmap_C, file = paste0("/import/home/share/mingxuan/fine_mapping/analysis/results/height_XMAP_K",K,"_out_", eur_data, "_", eas_data, "_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))

    # XMAP model without Omega background
    xmap_O0 <- XMAP::XMAP(simplify2array(list(R_eas, R_brit)), cbind(z_eas, z_eur), c(median(sumstat_EAS_i$N), median(sumstat_EUR_i$N)),
                          K = K, Omega = matrix(c(1e-30, 0, 0, 1e-30), 2, 2), Sig_E = c(c1, c2), tol = 1e-6,
                          maxIter = 200, estimate_residual_variance = F, estimate_prior_variance = T,
                          estimate_background_variance = F)
    saveRDS(xmap_O0, file = paste0("/import/home/share/mingxuan/fine_mapping/analysis/results/height_XMAP_Omega0_K",K,"_out_", eur_data, "_", eas_data, "_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))

    # XMAP model without C confounding adjustment
    xmap_C0 <- XMAP::XMAP(simplify2array(list(R_eas, R_brit)), cbind(z_eas, z_eur), c(median(sumstat_EAS_i$N), median(sumstat_EUR_i$N)),
                          K = K, Omega = OmegaHat, Sig_E = c(1, 1), tol = 1e-6,
                          maxIter = 200, estimate_residual_variance = F, estimate_prior_variance = T,
                          estimate_background_variance = F)
    saveRDS(xmap_C0, file = paste0("/import/home/share/mingxuan/fine_mapping/analysis/results/height_XMAP_C0_K",K,"_out_", eur_data, "_", eas_data, "_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))

    # XMAP model without C confounding adjustment
    xmap_0 <- XMAP::XMAP(simplify2array(list(R_eas, R_brit)), cbind(z_eas, z_eur), c(median(sumstat_EAS_i$N), median(sumstat_EUR_i$N)),
                         K = K, Omega = matrix(c(1e-30, 0, 0, 1e-30), 2, 2), Sig_E = c(1, 1), tol = 1e-6,
                         maxIter = 200, estimate_residual_variance = F, estimate_prior_variance = T,
                         estimate_background_variance = F)
    saveRDS(xmap_0, file = paste0("/import/home/share/mingxuan/fine_mapping/analysis/results/height_XMAP0_K",K,"_out_", eur_data, "_", eas_data, "_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))


    saveRDS(info, file = paste0("/import/home/share/mingxuan/fine_mapping/analysis/results/height_snpINFO_K",K,"_out_", eur_data, "_", eas_data, "_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))
    cat("CHR", chr, ", loci", i, "/", nrow(loci), " finished.\n")
  }
}

