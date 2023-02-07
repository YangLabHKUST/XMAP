library(data.table)
library(susieR)
source("/import/home/share/mingxuan/fine_mapping/utils.R")
# height_ukb <- fread("/home/share/mingxuan/Chinese_qc2/bolt/ukb_height_bolt.stats")
# height_wg <- fread("/home/share/mingxuan/Chinese_qc2/bolt/Chinese_height_bolt.stats")
#
# height_ukb$z <- height_ukb$BETA / height_ukb$SE
# height_wg$z <- height_wg$BETA / height_wg$SE

eur_data <- "UKB"
# eur_data <- "Sibship"

eas_data <- "WG"
# eas_data <- "BBJ"

thr_pip <- 0.8
K <- 10

snp_eas <- snp_eur <- snp_XMAP <- snp_XMAPC0 <- snp_XMAPO0 <- snp_XMAP0 <- data.frame()

out_eur <- out_eas <- out_xmap <- out_xmap0 <- out_xmapC0 <- out_xmapO0 <- data.frame()
for (chr in 2:22) {

  loci <- fread(paste0("/home/share/mingxuan/UKB_geno_finemap/chr", chr, "_eas_brit_afr.loci"))

  for (i in 1:nrow(loci)) {
    if (file.exists(paste0("/home/share/mingxuan/fine_mapping/analysis/results/height_snpINFO_K",K,"_out_", eur_data, "_", eas_data, "_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))) {

      snps_info <- readRDS(paste0("/home/share/mingxuan/fine_mapping/analysis/results/height_snpINFO_K",K,"_out_", eur_data, "_", eas_data, "_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))
      fit_XMAP <- readRDS(paste0("/home/share/mingxuan/fine_mapping/analysis/results/height_XMAP_K",K,"_out_", eur_data, "_", eas_data, "_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))
      fit_XMAPC0 <- readRDS(paste0("/home/share/mingxuan/fine_mapping/analysis/results/height_XMAP_C0_K",K,"_out_", eur_data, "_", eas_data, "_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))
      fit_eas <- readRDS(paste0("/home/share/mingxuan/fine_mapping/analysis/results/height_susie_K",K,"_out_", eas_data, "_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))
      fit_eur <- readRDS(paste0("/home/share/mingxuan/fine_mapping/analysis/results/height_susie_K",K,"_out_", eur_data, "_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))

      idx_1mb <- which(snps_info$position > loci$start[i] & snps_info$position < loci$end[i])
      snps_info <- snps_info[idx_1mb,]

      out_eur <- rbind(out_eur, data.frame(snps_info, PIP = fit_eur$pip[idx_1mb]))
      out_eas <- rbind(out_eas, data.frame(snps_info, PIP = fit_eas$pip[idx_1mb]))
      out_xmapC0 <- rbind(out_xmapC0, data.frame(snps_info, PIP = get_pip(fit_XMAPC0$gamma)[idx_1mb]))
      out_xmap <- rbind(out_xmap, data.frame(snps_info, PIP = get_pip(fit_XMAP$gamma)[idx_1mb]))

      pip_eas <- fit_eas$pip[idx_1mb]
      snp_eas <- rbind(snp_eas, data.frame(snps_info[pip_eas > thr_pip, 1:2], PIP = pip_eas[pip_eas > thr_pip]))

      pip_eur <- fit_eur$pip[idx_1mb]
      snp_eur <- rbind(snp_eur, data.frame(snps_info[pip_eur > thr_pip, 1:2], PIP = pip_eur[pip_eur > thr_pip]))



      pipC0 <- get_pip(fit_XMAPC0$gamma)[idx_1mb]
      snp_XMAPC0 <- rbind(snp_XMAPC0, data.frame(snps_info[pipC0 > thr_pip, 1:2], PIP = pipC0[pipC0 > thr_pip]))

      pip <- get_pip(fit_XMAP$gamma)[idx_1mb]
      snp_XMAP <- rbind(snp_XMAP, data.frame(snps_info[pip > thr_pip, 1:2], PIP = pip[pip > thr_pip]))

      # cat(i, "-th loci finished.\n")
    }
  }
  cat(chr, "-th chr finished.\n")
}

saveRDS(list(EAS = snp_eas, EUR = snp_eur, XMAP = snp_XMAP, XMAPC0 = snp_XMAPC0, XMAPO0 = snp_XMAPO0, XMAP0 = snp_XMAP0), file = paste0("/home/share/mingxuan/fine_mapping/analysis/results/height_snps_causal_K",K,"_thr", thr_pip, ".RDS"))
write.table(out_eur,file=paste0("/home/share/mingxuan/fine_mapping/analysis/results/height_",eur_data,"_K",K,"_PIP.txt"),col.names = T,row.names = F,sep="\t",quote=F)
write.table(out_eas,file=paste0("/home/share/mingxuan/fine_mapping/analysis/results/height_",eas_data,"_K",K,"_PIP.txt"),col.names = T,row.names = F,sep="\t",quote=F)

write.table(out_xmap,file=paste0("/home/share/mingxuan/fine_mapping/analysis/results/height_", eur_data, "_", eas_data, "_K",K,"_PIP.txt"),col.names = T,row.names = F,sep="\t",quote=F)
write.table(out_xmapC0,file=paste0("/home/share/mingxuan/fine_mapping/analysis/results/height_", eur_data, "_", eas_data, "_CI_K",K,"_PIP.txt"),col.names = T,row.names = F,sep="\t",quote=F)
