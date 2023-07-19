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
n1 <- 20000
n2 <- 20000
K_true <- 3
K <- 5
df <- 16

h1_omega <- 0.005
h2_omega <- 0.005
rho_omega <- 0 #0.8

enrich <- 50
rho_sigma <- 0

loci <- fread(paste0("/import/home/share/mingxuan/fine_mapping/loci_chr22_1mb.txt"))
i_loci <- 29
idx_include <- (501:1000) + loci$idx_start[i_loci] - 1
p <- length(idx_include)

bim_loci <- bim[bim$V4 >= loci$left[i_loci] & bim$V4 <= loci$right[i_loci],]
bp_start <- bim_loci$V4[idx_include[1]]
bp_end <- bim_loci$V4[idx_include[p]]
snps <- bim_loci$V2[idx_include]

# construct in-sample genotype for gwas and ld matrix in pop1
# n1_all <- c(5000, 10000, 15000, 20000)
# for (n1 in n1_all) {
#
#   ind_n1 <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_indiv_wg_in_sample_ref_n", n1, ".txt"))$x
#   iid_fid <- cbind(ind_n1, ind_n1)
#   write.table(iid_fid, paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_iid_fid_wg_in_sample_ref_n", n1, ".txt"), col.names = F, row.names = F, quote = F)
#
#   susiex_ld <- paste0("python SuSiEx/SuSiEx_LD.py --ref_file=/home/share/mingxuan/wegene_qc2/height_merge_qc2_chr22 --chr=22 --bp=", bp_start, ",", bp_end, " --plink=/import/home2/mcaiad/plink/plink --maf=0.0001")
#   susiex_ld <- paste0(susiex_ld, " --keep=/import/home/share/mingxuan/fine_mapping/simulation/data/simu_iid_fid_wg_in_sample_ref_n", n1, ".txt")
#   susiex_ld <- paste0(susiex_ld, " --ld_file=/home/share/mingxuan/fine_mapping/simulation/susiex/RunDirectory/LD1_n", n1)
#
#   system(susiex_ld)
#
#   system(paste0("mv /home/share/mingxuan/fine_mapping/simulation/susiex/RunDirectory/LD1_n",n1,".ld.gz /home/share/mingxuan/fine_mapping/simulation/susiex/RunDirectory/LD1_n",n1,".ld.gz.plink"))
#   R1 <- data.matrix(fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_R1_wg_in_sample_ref_n", n1, "_chr22_loci", i_loci, ".txt")))
#   R1 <- R1[idx_include, idx_include]
#   fwrite(R1,file=paste0("/home/share/mingxuan/fine_mapping/simulation/susiex/RunDirectory/LD1_n", n1,".ld.gz"),col.names = F,row.names = F,sep="\t")
#
# }
#
#
# n2_all <- c(10000, 15000, 20000)
# for (n2 in n2_all) {
#   ind_n2 <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_indiv_ukb_in_sample_ref_n", n2, ".txt"))$x
#   iid_fid <- cbind(ind_n2, ind_n2)
#   write.table(iid_fid, paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_iid_fid_ukb_in_sample_ref_n", n2, ".txt"), col.names = F, row.names = F, quote = F)
#
#   susiex_ld <- paste0("python SuSiEx/SuSiEx_LD.py --ref_file=/home/share/mingxuan/wegene_qc2/height_ukb_qc2_chr22 --chr=22 --bp=", bp_start, ",", bp_end, " --plink=/import/home2/mcaiad/plink/plink --maf=0.0001")
#   susiex_ld <- paste0(susiex_ld, " --keep=/import/home/share/mingxuan/fine_mapping/simulation/data/simu_iid_fid_ukb_in_sample_ref_n", n2, ".txt")
#   susiex_ld <- paste0(susiex_ld, " --ld_file=/home/share/mingxuan/fine_mapping/simulation/susiex/RunDirectory/LD2_n", n2)
#
#   system(susiex_ld)
#
#   system(paste0("mv /home/share/mingxuan/fine_mapping/simulation/susiex/RunDirectory/LD2_n",n2,".ld.gz /home/share/mingxuan/fine_mapping/simulation/susiex/RunDirectory/LD2_n",n2,".ld.gz.plink"))
#   R2 <- data.matrix(fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_R2_ukb_in_sample_ref_n", n2, "_chr22_loci", i_loci, ".txt")))
#   R2 <- R2[idx_include, idx_include]
#   fwrite(R2,file=paste0("/home/share/mingxuan/fine_mapping/simulation/susiex/RunDirectory/LD2_n", n2,".ld.gz"),col.names = F,row.names = F,sep="\t")
# }


for (i in 1:nrep) {
  z_scores1 <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_zs1_t", df, "_new_wholeCHR_loci_", i_loci, "_small_p", p, "_chr22_n", n1, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                            "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"))
  z_scores2 <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_zs2_t", df, "_new_wholeCHR_loci_", i_loci, "_small_p", p, "_chr22_n", n2, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                            "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"))
  ss1 <- data.frame(chr = 22, snp = snps, bp = bim_loci$V4[idx_include], a1 = bim_loci$V5[idx_include], a2 = bim_loci$V6[idx_include],
                    BETA = z_scores1$V2[loci$idx_left[i_loci]:loci$idx_right[i_loci]][idx_include] / sqrt(n1),
                    se = 1 / sqrt(n1), pval = 2 * pnorm(abs(z_scores1$V2[loci$idx_left[i_loci]:loci$idx_right[i_loci]][idx_include]), lower.tail = F))
  ss2 <- data.frame(chr = 22, snp = snps, bp = bim_loci$V4[idx_include], a1 = bim_loci$V5[idx_include], a2 = bim_loci$V6[idx_include],
                    BETA = z_scores2$V2[loci$idx_left[i_loci]:loci$idx_right[i_loci]][idx_include] / sqrt(n2),
                    se = 1 / sqrt(n2), pval = 2 * pnorm(abs(z_scores2$V2[loci$idx_left[i_loci]:loci$idx_right[i_loci]][idx_include]), lower.tail = F))

  write.table(ss1, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/susiex/RunDirectory/ss1_t", df, "_new_wholeCHR_loci_", i_loci, "_small_p", p, "_chr22_n", n1, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                 "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"), col.names = T, row.names = F, quote = F)

  write.table(ss2, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/susiex/RunDirectory/ss2_t", df, "_new_wholeCHR_loci_", i_loci, "_small_p", p, "_chr22_n", n2, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                 "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"), col.names = T, row.names = F, quote = F)

  susiex <- paste0("export NUMEXPR_NUM_THREADS=10; export OMP_NUM_THREADS=10; python SuSiEx/SuSiEx.py --pval_thresh=1 --min_purity=0 --max_iter=500 --out_dir=/import/home/share/mingxuan/fine_mapping/simulation/susiex/output --chr=22 --bp=", bp_start, ",", bp_end, " --chr_col=1,1 --snp_col=2,2 --bp_col=3,3 --a1_col=4,4 --a2_col=5,5 --eff_col=6,6 --se_col=7,7 --pval_col=8,8 --keep-ambig=True --mult-step=False --maf=0.0001 --full_out=True")
  susiex <- paste0(susiex, " --n_sig=", K)
  susiex <- paste0(susiex, " --sst_file=/import/home/share/mingxuan/fine_mapping/simulation/susiex/RunDirectory/ss1_t", df, "_new_wholeCHR_loci_", i_loci, "_small_p", p, "_chr22_n", n1, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                   "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt,/import/home/share/mingxuan/fine_mapping/simulation/susiex/RunDirectory/ss2_t", df, "_new_wholeCHR_loci_", i_loci, "_small_p", p, "_chr22_n", n2, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                   "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt")
  susiex <- paste0(susiex, " --ld_file=/import/home/share/mingxuan/fine_mapping/simulation/susiex/RunDirectory/LD1_n", n1, ",/import/home/share/mingxuan/fine_mapping/simulation/susiex/RunDirectory/LD2_n", n2)
  susiex <- paste0(susiex, " --n_gwas=", n1, ",", n2)
  susiex <- paste0(susiex, " --out_name=susiex_t", df, "_new_output_wholeCHR_loci_", i_loci, "_small_p", p, "_chr22_in_sample_n", n1, "_", n2,
                   "_Ktrue", K_true, "_K", K, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                   "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i)

  system(susiex)
  cat(i, "-th rep finished.\n")
}