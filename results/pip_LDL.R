library(data.table)
library(susieR)
source("/import/home/share/mingxuan/fine_mapping/utils.R")
thr_pip <- 0.99

snp_eas <- snp_afr <- snp_eur <- snp_XMAP <- snp_XMAPC0 <- snp_XMAPO0 <- snp_XMAP0 <- data.frame()

snp_eur_afr <- snp_eur_eas <- snp_afr_eas <- data.frame()

out_eur <- out_afr <- out_eas <- out_3pop <- out_3pop0 <- out_3popC0 <- out_3popO0 <- data.frame()
out_eur_afr <- out_eur_afr0 <- out_eur_afrC0 <- out_eur_afrO0 <- data.frame()
out_eur_eas <- out_eur_eas0 <- out_eur_easC0 <- out_eur_easO0 <- data.frame()
out_afr_eas <- out_afr_eas0 <- out_afr_easC0 <- out_afr_easO0 <- data.frame()
for (chr in c(1:22)) {

  loci <- fread(paste0("/home/share/mingxuan/UKB_geno_finemap/chr", chr, "_eas_brit_afr.loci"))

  for (i in 1:nrow(loci)) {
    if (file.exists(paste0("/home/share/mingxuan/fine_mapping/analysis/results/LDL_XMAP_Omega_out_ukb_afr_2pop_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))) {

      snps_info <- readRDS(paste0("/home/share/mingxuan/fine_mapping/analysis/results/LDL_snpINFO_out_ukb_afr_eas_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))
      fit_XMAP <- readRDS(paste0("/home/share/mingxuan/fine_mapping/analysis/results/LDL_XMAP_Omega_out_ukb_afr_eas_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))
      fit_XMAPC0 <- readRDS(paste0("/home/share/mingxuan/fine_mapping/analysis/results/LDL_XMAP_C0_out_ukb_afr_eas_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))
      fit_eur <- readRDS(paste0("/home/share/mingxuan/fine_mapping/analysis/results/LDL_susie_out_ukb_3pop_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))
      fit_afr <- readRDS(paste0("/home/share/mingxuan/fine_mapping/analysis/results/LDL_susie_out_afr_3pop_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))
      fit_eas <- readRDS(paste0("/home/share/mingxuan/fine_mapping/analysis/results/LDL_susie_out_eas_3pop_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))

      idx_1mb <- which(snps_info$position > loci$start[i] & snps_info$position < loci$end[i])
      snps_info <- snps_info[idx_1mb,]


      out_eur <- rbind(out_eur, data.frame(snps_info, PIP = fit_eur$pip[idx_1mb]))
      out_afr <- rbind(out_afr, data.frame(snps_info, PIP = fit_afr$pip[idx_1mb]))
      out_eas <- rbind(out_eas, data.frame(snps_info, PIP = fit_eas$pip[idx_1mb]))
      out_3popC0 <- rbind(out_3popC0, data.frame(snps_info, PIP = get_pip(fit_XMAPC0$gamma)[idx_1mb]))
      out_3pop <- rbind(out_3pop, data.frame(snps_info, PIP = get_pip(fit_XMAP$gamma)[idx_1mb]))



      # EUR_AFR
      fit_XMAP <- readRDS(paste0("/home/share/mingxuan/fine_mapping/analysis/results/LDL_XMAP_Omega_out_ukb_afr_2pop_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))
      fit_XMAPC0 <- readRDS(paste0("/home/share/mingxuan/fine_mapping/analysis/results/LDL_XMAP_C0_out_ukb_afr_2pop_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))

      out_eur_afrC0 <- rbind(out_eur_afrC0, data.frame(snps_info, PIP = get_pip(fit_XMAPC0$gamma)[idx_1mb]))
      out_eur_afr <- rbind(out_eur_afr, data.frame(snps_info, PIP = get_pip(fit_XMAP$gamma)[idx_1mb]))

      # EUR_EAS
      fit_XMAP <- readRDS(paste0("/home/share/mingxuan/fine_mapping/analysis/results/LDL_XMAP_Omega_out_ukb_eas_2pop_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))
      fit_XMAPC0 <- readRDS(paste0("/home/share/mingxuan/fine_mapping/analysis/results/LDL_XMAP_C0_out_ukb_eas_2pop_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))

      out_eur_easC0 <- rbind(out_eur_easC0, data.frame(snps_info, PIP = get_pip(fit_XMAPC0$gamma)[idx_1mb]))
      out_eur_eas <- rbind(out_eur_eas, data.frame(snps_info, PIP = get_pip(fit_XMAP$gamma)[idx_1mb]))

      # AFR_EAS
      fit_XMAP <- readRDS(paste0("/home/share/mingxuan/fine_mapping/analysis/results/LDL_XMAP_Omega_out_afr_eas_2pop_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))
      fit_XMAPC0 <- readRDS(paste0("/home/share/mingxuan/fine_mapping/analysis/results/LDL_XMAP_C0_out_afr_eas_2pop_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))

      out_afr_easC0 <- rbind(out_afr_easC0, data.frame(snps_info, PIP = get_pip(fit_XMAPC0$gamma)[idx_1mb]))
      out_afr_eas <- rbind(out_afr_eas, data.frame(snps_info, PIP = get_pip(fit_XMAP$gamma)[idx_1mb]))


      # cat(chr, "-th chr; ", i, "/", nrow(loci), " loci finished.\n")
    }
  }
  cat(chr, "-th chr finished.\n")
}

write.table(out_eur,file="/home/share/mingxuan/fine_mapping/analysis/results/LDL_EUR_PIP.txt",col.names = T,row.names = F,sep="\t",quote=F)
write.table(out_afr,file="/home/share/mingxuan/fine_mapping/analysis/results/LDL_AFR_PIP.txt",col.names = T,row.names = F,sep="\t",quote=F)
write.table(out_eas,file="/home/share/mingxuan/fine_mapping/analysis/results/LDL_EAS_PIP.txt",col.names = T,row.names = F,sep="\t",quote=F)

write.table(out_3pop,file="/home/share/mingxuan/fine_mapping/analysis/results/LDL_3pop_PIP.txt",col.names = T,row.names = F,sep="\t",quote=F)
write.table(out_3popC0,file="/home/share/mingxuan/fine_mapping/analysis/results/LDL_3pop_CI_PIP.txt",col.names = T,row.names = F,sep="\t",quote=F)

write.table(out_eur_afr,file="/home/share/mingxuan/fine_mapping/analysis/results/LDL_EUR_AFR_PIP.txt",col.names = T,row.names = F,sep="\t",quote=F)
write.table(out_eur_afrC0,file="/home/share/mingxuan/fine_mapping/analysis/results/LDL_EUR_AFR_CI_PIP.txt",col.names = T,row.names = F,sep="\t",quote=F)

write.table(out_eur_eas,file="/home/share/mingxuan/fine_mapping/analysis/results/LDL_EUR_EAS_PIP.txt",col.names = T,row.names = F,sep="\t",quote=F)
write.table(out_eur_easC0,file="/home/share/mingxuan/fine_mapping/analysis/results/LDL_EUR_EAS_CI_PIP.txt",col.names = T,row.names = F,sep="\t",quote=F)

write.table(out_afr_eas,file="/home/share/mingxuan/fine_mapping/analysis/results/LDL_AFR_EAS_PIP.txt",col.names = T,row.names = F,sep="\t",quote=F)
write.table(out_afr_easC0,file="/home/share/mingxuan/fine_mapping/analysis/results/LDL_AFR_EAS_CI_PIP.txt",col.names = T,row.names = F,sep="\t",quote=F)
