################################### Process fine mapping results ###################################
library(data.table)
trait <- "Hb"
out_eur <- fread(file = paste0("/home/share/mingxuan/fine_mapping/analysis/results/", trait, "_ukb_PIP.txt"), data.table = F)
out_eas <- fread(file = paste0("/home/share/mingxuan/fine_mapping/analysis/results/", trait, "_bbj_PIP.txt"), data.table = F)

out_xmap <- fread(file = paste0("/home/share/mingxuan/fine_mapping/analysis/results/", trait, "_ukb_bbj_PIP.txt"), data.table = F)

out_xmap <- data.frame(CHR = paste0("chr", out_xmap$chromosome), start = out_xmap$position, end = out_xmap$position + 1, PIP = out_xmap$PIP)

out_eur <- data.frame(CHR = paste0("chr", out_eur$chromosome), start = out_eur$position, end = out_eur$position + 1, PIP = out_eur$PIP)
out_eas <- data.frame(CHR = paste0("chr", out_eas$chromosome), start = out_eas$position, end = out_eas$position + 1, PIP = out_eas$PIP)


write.table(out_xmap, file = paste0("/home/share/mingxuan/fine_mapping/analysis/results/", trait, "_ukb_bbj.bed"), col.names = F, row.names = F, quote = F, sep = "\t")

write.table(out_eur, file = paste0("/home/share/mingxuan/fine_mapping/analysis/results/", trait, "_ukb.bed"), col.names = F, row.names = F, quote = F, sep = "\t")
write.table(out_eas, file = paste0("/home/share/mingxuan/fine_mapping/analysis/results/", trait, "_bbj.bed"), col.names = F, row.names = F, quote = F, sep = "\t")



library(SCAVENGE)
library(chromVAR)
library(gchromVAR)
# library(BuenColors)
library(SummarizedExperiment)
library(data.table)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
library(igraph)
library(ggplot2)

hemarda <- "/import/home/share/mingxuan/Data/scATAC-seq/scATAC-Healthy-Hematopoiesis-191120.rds"
SE_hema <- readRDS(hemarda)
SE_hema <- SE_hema[, -grep("Unk", SE_hema$BioClassification)]
SE_hema$BioClassification <- gsub(".*_", "", SE_hema$BioClassification)
SE_hema <- addGCBias(SE_hema, genome = BSgenome.Hsapiens.UCSC.hg19)
SE_hema_bg <- getBackgroundPeaks(SE_hema, niterations = 200)

# peak_by_cell_mat <- assay(SE_hema)
# tfidf_mat <- tfidf(bmat=peak_by_cell_mat, mat_binary=TRUE, TF=TRUE, log_TF=TRUE)
# lsi_mat <- do_lsi(tfidf_mat, dims=30)
# mutualknn30 <- getmutualknn(lsi_mat, 30)
# saveRDS(mutualknn30,file="/import/home/share/mingxuan/Data/scATAC-seq/scATAC-Healthy-Hematopoiesis-191120_filterUnk_mknn30.RDS")
mutualknn30 <- readRDS("/import/home/share/mingxuan/Data/scATAC-seq/scATAC-Healthy-Hematopoiesis-191120_filterUnk_mknn30.RDS")


# trait <- "Baso"

bbj_file <- paste0("/home/share/mingxuan/fine_mapping/analysis/results/", trait, "_bbj.bed")
ukb_file <- paste0("/home/share/mingxuan/fine_mapping/analysis/results/", trait, "_ukb.bed")
xmap_file <- paste0("/home/share/mingxuan/fine_mapping/analysis/results/", trait, "_ukb_bbj.bed")

bbj_import <- importBedScore(rowRanges(SE_hema), bbj_file, colidx = 4)
ukb_import <- importBedScore(rowRanges(SE_hema), ukb_file, colidx = 4)
xmap_import <- importBedScore(rowRanges(SE_hema), xmap_file, colidx = 4)

SE_hema_bbj <- computeWeightedDeviations(SE_hema, bbj_import, background_peaks = SE_hema_bg)
SE_hema_ukb <- computeWeightedDeviations(SE_hema, ukb_import, background_peaks = SE_hema_bg)
SE_hema_xmap <- computeWeightedDeviations(SE_hema, xmap_import, background_peaks = SE_hema_bg)

z_score_bbj <- data.frame(colData(SE_hema), z_score = t(assays(SE_hema_bbj)[["z"]]) %>% c)
z_score_ukb <- data.frame(colData(SE_hema), z_score = t(assays(SE_hema_ukb)[["z"]]) %>% c)
z_score_xmap <- data.frame(colData(SE_hema), z_score = t(assays(SE_hema_xmap)[["z"]]) %>% c)
# head(z_score_bbj)
# head(z_score_ukb)
# head(z_score_xmap)


for (fm_method in c("bbj", "ukb", "xmap")) {
  z_score_trait <- get(paste0("z_score_", fm_method))

  seed_idx <- seedindex(z_score_trait$z_score, 0.05)
  scale_factor <- cal_scalefactor(z_score = z_score_trait$z_score, 0.01)

  np_score <- randomWalk_sparse(intM = mutualknn30, rownames(mutualknn30)[seed_idx], gamma = 0.05)

  omit_idx <- np_score == 0
  sum(omit_idx)

  np_score <- np_score[!omit_idx]
  TRS <- np_score %>%
    capOutlierQuantile(., 0.95) %>%
    max_min_scale
  TRS <- TRS * scale_factor
  trait_mat <- data.frame(z_score_trait[!omit_idx,], seed_idx[!omit_idx], np_score, TRS)

  trait_permu <- get_sigcell_simple(knn_sparse_mat = mutualknn30[!omit_idx, !omit_idx], seed_idx = trait_mat$seed_idx, topseed_npscore = trait_mat$np_score, permutation_times = 1000, true_cell_significance = 0.05, rda_output = F, mycores = 8, rw_gamma = 0.05)
  trait_mat <- data.frame(trait_mat, trait_permu)
  write.table(trait_mat, file = paste0("/home/share/mingxuan/fine_mapping/analysis/results/scavenge_output_", trait, "_", fm_method, ".txt"), col.names = T, row.names = F, sep = "\t", quote = F)
}
