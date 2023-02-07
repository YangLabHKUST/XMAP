library(susieR)
library(data.table)
library(Matrix)
library(XMAP)
library(RhpcBLASctl)
set.seed(1)
blas_set_num_threads(20)



# Step 1: estimate confounding and polygenic parameters
sumstat_GLGC <- fread("/import/home/share/sumstats/format/LDL_EUR_GLGC_summary_format.txt")
sumstat_EUR <- fread("/import/home/share/sumstats/format/LDL_EUR_GLGC_summary_format.txt")
sumstat_EUR$Z <- sumstat_EUR$beta/sumstat_EUR$se
# sumstat_EUR <- fread("/import/home/share/sumstats/format/LDL_allSNPs_UKBNealLab_summary_format.txt")
sumstat_AFR <- fread("/import/home/share/sumstats/format/LDL_AFR_GLGC_summary_format.txt")
sumstat_EAS <- fread("/import/home/share/sumstats/format/LDL_EAS_GLGC_summary_format.txt")

ldscore <- data.frame()
for (chr in 1:22) {
  ldscore_chr <- fread(paste0("/import/home/share/mingxuan/UKB_geno_finemap/LDscore_eas_brit_afr_chr", chr, ".txt"))
  ldscore <- rbind(ldscore, ldscore_chr)
  cat("CHR", chr, "\n")
}

# remove ambiguous SNPs
idx_amb <- which(ldscore$allele1 == comple(ldscore$allele2))
ldscore <- ldscore[-idx_amb,]


# overlap SNPs
snps <- Reduce(intersect, list(ldscore$rsid, sumstat_AFR$SNP, sumstat_EUR$SNP, sumstat_EAS$SNP))
sumstat_AFR <- sumstat_AFR[match(snps, sumstat_AFR$SNP),]
sumstat_EUR <- sumstat_EUR[match(snps, sumstat_EUR$SNP),]
sumstat_EAS <- sumstat_EAS[match(snps, sumstat_EAS$SNP),]
ldscore <- ldscore[match(snps, ldscore$rsid),]


# flip alleles
z_afr <- sumstat_AFR$beta / sumstat_AFR$se
z_eas <- sumstat_EAS$beta / sumstat_EAS$se
z_eur <- sumstat_EUR$Z

idx_flip <- which(sumstat_AFR$A1 != ldscore$allele1 & sumstat_AFR$A1 != comple(ldscore$allele1))
z_afr[idx_flip] <- -z_afr[idx_flip]

idx_flip <- which(sumstat_EUR$A1 != ldscore$allele1 & sumstat_EUR$A1 != comple(ldscore$allele1))
z_eur[idx_flip] <- -z_eur[idx_flip]

idx_flip <- which(sumstat_EAS$A1 != ldscore$allele1 & sumstat_EAS$A1 != comple(ldscore$allele1))
z_eas[idx_flip] <- -z_eas[idx_flip]


idx1 <- which(z_afr^2 < 30 & z_eur^2 < 30 & z_eas^2 < 30)
ld_eas_w <- 1 / sapply(ldscore$EAS, function(x) max(x, 1))
ld_afr_w <- 1 / sapply(ldscore$AFR, function(x) max(x, 1))
ld_eur_w <- 1 / sapply(ldscore$EUR, function(x) max(x, 1))

# AFR-EUR
fit_step1 <- estimate_gc(data.frame(Z = z_afr[idx1], N = sumstat_AFR$N[idx1]), data.frame(Z = z_eur[idx1], N = sumstat_EUR$N[idx1]),
                         ldscore$AFR[idx1], ldscore$EUR[idx1], ldscore$AFR_EUR[idx1],
                         reg_w1 = ld_afr_w[idx1], reg_w2 = ld_eur_w[idx1], reg_wx = sqrt(ld_afr_w[idx1] * ld_eur_w[idx1]),
                         constrain_intercept = F)
fit_step2 <- estimate_gc(data.frame(Z = z_afr, N = sumstat_AFR$N), data.frame(Z = z_eur, N = sumstat_EUR$N),
                         ldscore$AFR, ldscore$EUR, ldscore$AFR_EUR,
                         reg_w1 = ld_afr_w, reg_w2 = ld_eur_w, reg_wx = sqrt(ld_afr_w * ld_eur_w),
                         constrain_intercept = T, fit_step1$tau1$coefs[1], fit_step1$tau2$coefs[1], fit_step1$theta$coefs[1])

# $tau1
# [1] 1.066208e+00 3.882615e-08
# $tau2
# [1] 1.094689e+00 1.931002e-08
# $theta
# [1] 2.235382e-02 1.975517e-08

# EAS-EUR
fit_step1 <- estimate_gc(data.frame(Z = z_eas[idx1], N = sumstat_EAS$N[idx1]), data.frame(Z = z_eur[idx1], N = sumstat_EUR$N[idx1]),
                         ldscore$EAS[idx1], ldscore$EUR[idx1], ldscore$EUR_EAS[idx1],
                         reg_w1 = ld_eas_w[idx1], reg_w2 = ld_eur_w[idx1], reg_wx = sqrt(ld_eas_w[idx1] * ld_eur_w[idx1]),
                         constrain_intercept = F)
fit_step2 <- estimate_gc(data.frame(Z = z_eas, N = sumstat_EAS$N), data.frame(Z = z_eur, N = sumstat_EUR$N),
                         ldscore$EAS, ldscore$EUR, ldscore$EUR_EAS,
                         reg_w1 = ld_eas_w, reg_w2 = ld_eur_w, reg_wx = sqrt(ld_eas_w * ld_eur_w),
                         constrain_intercept = T, fit_step1$tau1$coefs[1], fit_step1$tau2$coefs[1], fit_step1$theta$coefs[1])

# $tau1
# [1] 1.036577e+00 2.370618e-08
# $tau2
# [1] 1.094689e+00 1.931002e-08
# $theta
# [1] 1.754128e-02 1.359512e-08


# EAS-AFR
fit_step1 <- estimate_gc(data.frame(Z = z_eas[idx1], N = sumstat_EAS$N[idx1]), data.frame(Z = z_afr[idx1], N = sumstat_AFR$N[idx1]),
                         ldscore$EAS[idx1], ldscore$AFR[idx1], ldscore$AFR_EAS[idx1],
                         reg_w1 = ld_eas_w[idx1], reg_w2 = ld_afr_w[idx1], reg_wx = sqrt(ld_eas_w[idx1] * ld_afr_w[idx1]),
                         constrain_intercept = F)
fit_step2 <- estimate_gc(data.frame(Z = z_eas, N = sumstat_EAS$N), data.frame(Z = z_afr, N = sumstat_AFR$N),
                         ldscore$EAS, ldscore$AFR, ldscore$AFR_EAS,
                         reg_w1 = ld_eas_w, reg_w2 = ld_afr_w, reg_wx = sqrt(ld_eas_w * ld_afr_w),
                         constrain_intercept = T, fit_step1$tau1$coefs[1], fit_step1$tau2$coefs[1], fit_step1$theta$coefs[1])

# $tau1
# [1] 1.036577e+00 2.370618e-08
# $tau2
# [1] 1.066208e+00 3.882615e-08
# $theta
# [1] 6.959741e-03 2.102276e-08





# Step 2: Evaluate SNP PIP given confounding and polygenic parameters
library(susieR)
library(data.table)
library(Matrix)
library(XMAP)
library(RhpcBLASctl)
set.seed(1)
blas_set_num_threads(15)
source("/import/home/share/mingxuan/fine_mapping/utils.R")
source("/import/home/share/mingxuan/fine_mapping/XMAP.R")
source("/home/share/mingxuan/taam/code/sldxr.R")

# sumstat_EUR <- fread("/home/share/sumstats/format/LDL_EUR_GLGC_summary_format.txt")
sumstat_EUR <- fread("/home/share/sumstats/format/LDL_allSNPs_UKBNealLab_summary_format.txt")
sumstat_AFR <- fread("/home/share/sumstats/format/LDL_AFR_GLGC_summary_format.txt")
sumstat_EAS <- fread("/home/share/sumstats/format/LDL_EAS_GLGC_summary_format.txt")
OmegaHat <- diag(c(1.931002e-08, 3.882615e-08, 2.370618e-08)) # EUR AFR EAS
OmegaHat[1, 2] <- 1.975517e-08
OmegaHat[1, 3] <- 1.359512e-08
OmegaHat[2, 3] <- 2.102276e-08
OmegaHat[lower.tri(OmegaHat)] <- OmegaHat[upper.tri(OmegaHat)]

c1 <- 1.094689
c2 <- 1.066208
c3 <- 1.036577

chr_all <- 12
for (chr in chr_all) {


  loci <- fread(paste0("/home/share/mingxuan/UKB_geno_finemap/chr", chr, "_eas_brit_afr.loci"))
  for (i in 1:nrow(loci)) {

    info <- fread(paste0("/home/share/mingxuan/UKB_geno_finemap/chr", chr, "_", loci$left[i], "_", loci$right[i], "_eas.info"))
    info_eur_afr <- fread(paste0("/home/share/mingxuan/UKB_geno_finemap/chr", chr, "_", loci$left[i], "_", loci$right[i], ".info"))

    # remove ambiguous SNPs
    idx_amb <- which(info$allele1 == comple(info$allele2))
    info <- info[-idx_amb,]


    snps <- Reduce(intersect, list(info$rsid, sumstat_AFR$SNP, sumstat_EUR$SNP, sumstat_EAS$SNP))

    if (length(snps) < 10) next
    # overlap SNPs
    sumstat_EAS_i <- sumstat_EAS[match(snps, sumstat_EAS$SNP),]
    sumstat_AFR_i <- sumstat_AFR[match(snps, sumstat_AFR$SNP),]
    sumstat_EUR_i <- sumstat_EUR[match(snps, sumstat_EUR$SNP),]
    if (min(c(2 * pnorm(abs(sumstat_AFR_i$beta / sumstat_AFR_i$se), lower.tail = F),
              2 * pnorm(abs(sumstat_EUR_i$Z), lower.tail = F),
              2 * pnorm(abs(sumstat_EAS_i$beta / sumstat_EAS_i$se), lower.tail = F))) > 5e-8) {
      cat(i, "-th loci has no significant snps.\n")
      next
    }
    R_eas <- readMM(paste0("/home/share/mingxuan/UKB_geno_finemap/chr", chr, "_", loci$left[i], "_", loci$right[i], "_eas.mtx.gz"))
    R_eas <- as.matrix(R_eas + t(R_eas))
    R_eas <- R_eas[-idx_amb, -idx_amb]
    idx_eas <- match(snps, info$rsid)
    R_eas <- R_eas[idx_eas, idx_eas]
    info <- info[idx_eas,]

    R_afr <- readMM(paste0("/home/share/mingxuan/UKB_geno_finemap/chr", chr, "_", loci$left[i], "_", loci$right[i], "_afr.mtx.gz"))
    R_afr <- as.matrix(R_afr + t(R_afr))
    idx_afr <- match(snps, info_eur_afr$rsid)
    R_afr <- R_afr[idx_afr, idx_afr]

    R_brit <- readMM(paste0("/home/share/mingxuan/UKB_geno_finemap/chr", chr, "_", loci$left[i], "_", loci$right[i], "_brit.mtx.gz"))
    R_brit <- as.matrix(R_brit + t(R_brit))
    idx_brit <- match(snps, info_eur_afr$rsid)
    R_brit <- R_brit[idx_brit, idx_brit]


    # remove sample size outliers

    idx_outlier_EUR <- which(sumstat_EUR_i$N < 0.7 * median(sumstat_EUR_i$N))

    idx_outlier_AFR <- which(sumstat_AFR_i$N < 0.7 * median(sumstat_AFR_i$N))

    idx_outlier_EAS <- which(sumstat_EAS_i$N < 0.7 * median(sumstat_EAS_i$N))


    idx_outlier <- unique(c(idx_outlier_EUR, idx_outlier_AFR, idx_outlier_EAS))
    if (length(idx_outlier) > 0) {

      snps <- snps[-idx_outlier]
      sumstat_AFR_i <- sumstat_AFR_i[-idx_outlier,]
      sumstat_EUR_i <- sumstat_EUR_i[-idx_outlier,]
      sumstat_EAS_i <- sumstat_EAS_i[-idx_outlier,]
      info <- info[-idx_outlier,]
      R_eas <- R_eas[-idx_outlier, -idx_outlier]
      R_afr <- R_afr[-idx_outlier, -idx_outlier]
      R_brit <- R_brit[-idx_outlier, -idx_outlier]
    }

    # flip alleles
    z_afr <- sumstat_AFR_i$beta / sumstat_AFR_i$se
    z_eas <- sumstat_EAS_i$beta / sumstat_EAS_i$se
    z_eur <- sumstat_EUR_i$Z

    idx_flip <- which(sumstat_EAS_i$A1 != info$allele1 & sumstat_EAS_i$A1 != comple(info$allele1))
    z_eas[idx_flip] <- -z_eas[idx_flip]

    idx_flip <- which(sumstat_AFR_i$A1 != info$allele1 & sumstat_AFR_i$A1 != comple(info$allele1))
    z_afr[idx_flip] <- -z_afr[idx_flip]

    idx_flip <- which(sumstat_EUR_i$A1 != info$allele1 & sumstat_EUR_i$A1 != comple(info$allele1))
    z_eur[idx_flip] <- -z_eur[idx_flip]

    susie_eas <- susie_rss(z_eas, R_eas,
                           L = 10,
                           verbose = TRUE,
                           coverage = 0.9,
                           min_abs_corr = 0.1,
                           check_prior = F)
    saveRDS(susie_eas, file = paste0("/home/share/mingxuan/fine_mapping/analysis/results/LDL_susie_out_eas_3pop_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))
    susie_plot(susie_eas, y = "PIP", b = (info$rsid == "rs900776"), main = "GWAS-GLGC_EAS")


    susie_afr <- susie_rss(z_afr, R_afr,
                           L = 10,
                           verbose = TRUE,
                           coverage = 0.9,
                           min_abs_corr = 0.1,
                           check_prior = F)
    saveRDS(susie_afr, file = paste0("/home/share/mingxuan/fine_mapping/analysis/results/LDL_susie_out_afr_3pop_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))
    susie_plot(susie_afr, y = "PIP", b = (info$rsid == "rs900776"), main = "GWAS-GLGC_AFR")


    susie_eur <- susie_rss(z_eur, R_brit,
                           L = 10,
                           verbose = TRUE,
                           coverage = 0.9,
                           min_abs_corr = 0.1,
                           check_prior = F)
    saveRDS(susie_eur, file = paste0("/home/share/mingxuan/fine_mapping/analysis/results/LDL_susie_out_ukb_3pop_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))
    susie_plot(susie_eur, y = "PIP", b = (info$rsid == "rs900776"), main = "GWAS-UKB_EUR")


    xmap_C <- XMAP::XMAP(simplify2array(list(R_brit, R_afr, R_eas)), cbind(z_eur, z_afr, z_eas), c(median(sumstat_EUR_i$N), median(sumstat_AFR_i$N), median(sumstat_EAS_i$N)),
                         K = 10, Omega = OmegaHat, Sig_E = c(c1, c2, c3), tol = 1e-6,
                         maxIter = 200, estimate_residual_variance = T, estimate_prior_variance = T,
                         estimate_background_variance = F)
    saveRDS(xmap_C, file = paste0("/home/share/mingxuan/fine_mapping/analysis/results/LDL_XMAP_Omega_out_ukb_afr_eas_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))

    # XMAP model without C confounding adjustment
    xmap_C0 <- XMAP::XMAP(simplify2array(list(R_brit, R_afr, R_eas)), cbind(z_eur, z_afr, z_eas), c(median(sumstat_EUR_i$N), median(sumstat_AFR_i$N), median(sumstat_EAS_i$N)),
                          K = 10, Omega = OmegaHat, Sig_E = c(1, 1, 1), tol = 1e-6,
                          maxIter = 200, estimate_residual_variance = F, estimate_prior_variance = T,
                          estimate_background_variance = F)
    saveRDS(xmap_C0, file = paste0("/home/share/mingxuan/fine_mapping/analysis/results/LDL_XMAP_C0_out_ukb_afr_eas_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))

    saveRDS(info, file = paste0("/home/share/mingxuan/fine_mapping/analysis/results/LDL_snpINFO_out_ukb_afr_eas_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))

    cat("CHR", chr, ", loci", i, "/", nrow(loci), " finished.\n")
  }
}

