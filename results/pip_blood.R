library(data.table)
library(susieR)
source("/import/home/share/mingxuan/fine_mapping/utils.R")
trait <- "MCV"


thr_pip <- 0.99

snp_eas <- snp_eur <- snp_XMAP <- data.frame()

out_eur <- out_eas <- out_xmap <- data.frame()
for (chr in 1:22) {

  loci <- fread(paste0("/home/share/mingxuan/UKB_geno_finemap/chr", chr, "_eas_brit_afr.loci"))

  for (i in 1:nrow(loci)) {
    if (file.exists(paste0("/home/share/mingxuan/fine_mapping/analysis/results/",trait,"_snpINFO_out_ukb_eas_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))) {

      snps_info <- readRDS(paste0("/home/share/mingxuan/fine_mapping/analysis/results/",trait,"_snpINFO_out_ukb_eas_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))
      fit_XMAP <- readRDS(paste0("/home/share/mingxuan/fine_mapping/analysis/results/",trait,"_XMAP_Omega_out_ukb_eas_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))
      fit_eas <- readRDS(paste0("/home/share/mingxuan/fine_mapping/analysis/results/",trait,"_susie_out_eas_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))
      fit_eur <- readRDS(paste0("/home/share/mingxuan/fine_mapping/analysis/results/",trait,"_susie_out_ukb_chr", chr, "_", loci$left[i], "_", loci$right[i], ".RDS"))

      idx_1mb <- which(snps_info$position > loci$start[i] & snps_info$position < loci$end[i])
      snps_info <- snps_info[idx_1mb,]

      out_eur <- rbind(out_eur, data.frame(snps_info, PIP = fit_eur$pip[idx_1mb]))
      out_eas <- rbind(out_eas, data.frame(snps_info, PIP = fit_eas$pip[idx_1mb]))
      out_xmap <- rbind(out_xmap, data.frame(snps_info, PIP = get_pip(fit_XMAP$gamma)[idx_1mb]))


      # cat(i, "-th loci finished.\n")
    }
  }
  cat(chr, "-th chr finished.\n")
}
write.table(out_eur,file=paste0("/home/share/mingxuan/fine_mapping/analysis/results/",trait,"_ukb_PIP.txt"),col.names = T,row.names = F,sep="\t",quote=F)
write.table(out_eas,file=paste0("/home/share/mingxuan/fine_mapping/analysis/results/",trait,"_bbj_PIP.txt"),col.names = T,row.names = F,sep="\t",quote=F)

write.table(out_xmap,file=paste0("/home/share/mingxuan/fine_mapping/analysis/results/",trait,"_ukb_bbj_PIP.txt"),col.names = T,row.names = F,sep="\t",quote=F)



