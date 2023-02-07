#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

library(data.table)
library(susieR)
library(XMAP)

trait <- args[1]#"Hb"

sumstat_EUR <- fread(paste0("/home/share/sumstats/format/",trait,"_allSNPs_UKBNealLab_summary_format.txt"))
sumstat_EAS <- fread(paste0("/home/share/sumstats/format/",trait,"_BBJ_allSNPs_summary_format.txt"))


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
snps <- Reduce(intersect, list(ldscore$rsid, sumstat_EAS$SNP, sumstat_EUR$SNP))
sumstat_EAS <- sumstat_EAS[match(snps, sumstat_EAS$SNP),]
sumstat_EUR <- sumstat_EUR[match(snps, sumstat_EUR$SNP),]
ldscore <- ldscore[match(snps, ldscore$rsid),]


# flip alleles
z_eas <- sumstat_EAS$Z
# z_eur <- sumstat_EUR$beta / sumstat_EUR$se
z_eur <- sumstat_EUR$Z

idx_flip <- which(sumstat_EAS$A1 != ldscore$allele1 & sumstat_EAS$A1 != comple(ldscore$allele1))
z_eas[idx_flip] <- -z_eas[idx_flip]

idx_flip <- which(sumstat_EUR$A1 != ldscore$allele1 & sumstat_EUR$A1 != comple(ldscore$allele1))
z_eur[idx_flip] <- -z_eur[idx_flip]


idx1 <- which(z_eas^2 < 30 & z_eur^2 < 30)

ld_eas_w <- 1 / sapply(ldscore$EAS, function(x) max(x, 1))
ld_eur_w <- 1 / sapply(ldscore$EUR, function(x) max(x, 1))
fit_step1 <- estimate_gc(data.frame(Z = z_eas[idx1], N = sumstat_EAS$N[idx1]), data.frame(Z = z_eur[idx1], N = sumstat_EUR$N[idx1]),
                         ldscore$EAS[idx1], ldscore$EUR[idx1], ldscore$EUR_EAS[idx1],
                         reg_w1 = ld_eas_w[idx1], reg_w2 = ld_eur_w[idx1], reg_wx = sqrt(ld_eas_w[idx1] * ld_eur_w[idx1]),
                         constrain_intercept = F)


fit_step2 <- estimate_gc(data.frame(Z = z_eas, N = sumstat_EAS$N), data.frame(Z = z_eur, N = sumstat_EUR$N),
                         ldscore$EAS, ldscore$EUR, ldscore$EUR_EAS,
                         reg_w1 = ld_eas_w, reg_w2 = ld_eur_w, reg_wx = sqrt(ld_eas_w * ld_eur_w),
                         constrain_intercept = T, fit_step1$tau1$coefs[1], fit_step1$tau2$coefs[1], fit_step1$theta$coefs[1])
sink(paste0("/home/share/mingxuan/fine_mapping/analysis/results/",trait,"_ldsc_out.txt"))
print(fit_step1)
print(fit_step2)
sink()
################################################
# $tau1
# [1] 1.120054e+00 3.022511e-08
# $tau2
# [1] 1.270350e+00 4.617689e-08
# $theta
# [1] 6.611797e-02 3.055873e-08
################################################

library(data.table)
library(RhpcBLASctl)
library(susieR)
library(Matrix)
blas_set_num_threads(20)
# library(susieR)
source("/import/home/share/mingxuan/fine_mapping/utils.R")
source("/import/home/share/mingxuan/fine_mapping/XMAP.R")
source("/home/share/mingxuan/taam/code/sldxr.R")
sumstat_EUR <- fread(paste0("/home/share/sumstats/format/",trait,"_allSNPs_UKBNealLab_summary_format.txt"))
sumstat_EAS <- fread(paste0("/home/share/sumstats/format/",trait,"_BBJ_allSNPs_summary_format.txt"))

c1 <- fit_step2$tau1[1]
c2 <- fit_step2$tau2[1]
Omega1 <- fit_step2$tau1[2]
Omega2 <- fit_step2$tau2[2]
Omega12 <- fit_step2$theta[2]
OmegaHat <- OmegaHat0 <- diag(c(Omega1, Omega2))
OmegaHat[1, 2] <- OmegaHat[2, 1] <- Omega12
###################################################################


chr_all <- 1:22 #as.numeric(args[1])
for (chr in chr_all) {
  loci <- fread(paste0("/home/share/mingxuan/UKB_geno_finemap/chr", chr, "_eas_brit_afr.loci"))
  for (i in 1:nrow(loci)) {
    # for (i in 81) {

    info <- fread(paste0("/home/share/mingxuan/UKB_geno_finemap/chr", chr, "_", loci$left[i], "_", loci$right[i], "_eas.info"))
    info_eur_afr <- fread(paste0("/home/share/mingxuan/UKB_geno_finemap/chr", chr, "_", loci$left[i], "_", loci$right[i], ".info"))
    # remove ambiguous SNPs
    idx_amb <- which(info$allele1 == comple(info$allele2))
    info <- info[-idx_amb,]


    snps <- Reduce(intersect, list(info$rsid, sumstat_EAS$SNP, sumstat_EUR$SNP))

    if (length(snps) < 10) next
    # overlap SNPs
    sumstat_EAS_i <- sumstat_EAS[match(snps, sumstat_EAS$SNP),]
    sumstat_EUR_i <- sumstat_EUR[match(snps, sumstat_EUR$SNP),]
    if (min(c(2 * pnorm(abs(sumstat_EUR_i$Z), lower.tail = F),
              2 * pnorm(abs(sumstat_EAS_i$Z), lower.tail = F))) > 5e-8) {
      cat(i, "-th loci has no significant snps.\n")
      next
    }
    R_eas <- readMM(paste0("/home/share/mingxuan/UKB_geno_finemap/chr", chr, "_", loci$left[i], "_", loci$right[i], "_eas.mtx.gz"))
    R_eas <- as.matrix(R_eas + t(R_eas))
    R_eas <- R_eas[-idx_amb, -idx_amb]
    idx_eas <- match(snps, info$rsid)
    R_eas <- R_eas[idx_eas, idx_eas]
    info <- info[idx_eas,]

    R_brit <- readMM(paste0("/home/share/mingxuan/UKB_geno_finemap/chr", chr, "_", loci$left[i], "_", loci$right[i], "_brit.mtx.gz"))
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
                           L = 10,
                           verbose = TRUE,
                           coverage = 0.9,
                           min_abs_corr = 0.1,
                           check_prior = F)
    saveRDS(susie_eas, file = paste0("/home/share/mingxuan/fine_mapping/analysis/results/",trait,"_susie_out_eas_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))
    susie_plot(susie_eas, y = "PIP",  main = "SuSiE-WeGene_EAS")
    susie_plot(z_eas, y = "z",  main = "SuSiE-WeGene_EAS")

    susie_eur <- susie_rss(z_eur, R_brit,
                           L = 10,
                           verbose = TRUE,
                           coverage = 0.9,
                           min_abs_corr = 0.1,
                           check_prior = F)
    saveRDS(susie_eur, file = paste0("/home/share/mingxuan/fine_mapping/analysis/results/",trait,"_susie_out_ukb_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))
    susie_plot(susie_eur, y = "PIP",  main = "SuSiE-UKB_EUR")
    susie_plot(z_eur, y = "z", main = "GWAS-UKB_EUR")



    # full XMAP model
    xmap_C <- XMAP::XMAP(simplify2array(list(R_eas, R_brit)), cbind(z_eas, z_eur), c(median(sumstat_EAS_i$N), median(sumstat_EUR_i$N)),
                         K = 10, Omega = OmegaHat, Sig_E = c(c1, c2), tol = 1e-6,
                         maxIter = 200, estimate_residual_variance = F, estimate_prior_variance = T,
                         estimate_background_variance = F)
    saveRDS(xmap_C, file = paste0("/home/share/mingxuan/fine_mapping/analysis/results/",trait,"_XMAP_Omega_out_ukb_eas_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))



    saveRDS(info, file = paste0("/home/share/mingxuan/fine_mapping/analysis/results/",trait,"_snpINFO_out_ukb_eas_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))
    cat("CHR", chr, ", loci", i, "/", nrow(loci), " finished.\n")
  }
}

