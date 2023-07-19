#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
########################################################### pop1 ###########################################################
library(VCM)
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
n1 <- 10000#as.numeric(args[1]) #10000
K_true <- 1#as.numeric(args[2]) #1
K <- 5

h1_omega <- 0.005
h2_omega <- 0.005
rho_omega <- 0.8

enrich <- 50
rho_sigma <- 0

pc_var1 <- 0.05
pc_var2 <- 0.2

loci <- fread(paste0("/home/share/mingxuan/fine_mapping/loci_chr22_1mb.txt"))
i_loci <- 29
idx_include <- (501:1000) + loci$idx_start[i_loci] - 1
p <- length(idx_include)

bim_loci <- bim[bim$V4 >= loci$left[i_loci] & bim$V4 <= loci$right[i_loci],]
bp_start <- bim_loci$V4[idx_include[1]]
bp_end <- bim_loci$V4[idx_include[p]]
snps <- bim_loci$V2[idx_include]

R1 <- data.matrix(fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_R1_wg_in_sample_ref_n", n1, "_chr22_loci", i_loci, ".txt")))
R1 <- R1[idx_include, idx_include]
fwrite(R1, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/susieinf/RunDirectory/R1_ref_n", n1, "_chr22_loci", i_loci, ".ld.gz"), col.names = F, quote = F, sep = "\t")
Beta1 <- data.matrix(fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_loci_", i_loci, "_small_p", p, "_Beta1_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                  "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".txt")))


for (i in 1:nrep) {
  z_scores1 <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_zs1_inflate_pc1_", pc_var1, "_corrected_loci_", i_loci, "_small_p", p, "_chr22_n", n1, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                            "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"))

  zs1 <- z_scores1$V2[loci$idx_left[i_loci]:loci$idx_right[i_loci]][idx_include]

  fwrite(data.frame(Z = zs1), file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/susieinf/RunDirectory/simu_zs1_inflate_susieinf_loci_", i_loci, "_small_p", p, "_chr22_n", n1, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                            "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"))

  susieinf1 <- "python fine-mapping-inf/run_fine_mapping.py --z-col-name Z --save-tsv"
  susieinf1 <- paste0(susieinf1, " --sumstats /import/home/share/mingxuan/fine_mapping/simulation/susieinf/RunDirectory/simu_zs1_inflate_susieinf_loci_", i_loci, "_small_p", p, "_chr22_n", n1, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                      "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt")
  susieinf1 <- paste0(susieinf1, " --num-sparse-effects ", K)
  susieinf1 <- paste0(susieinf1, " --ld-file /import/home/share/mingxuan/fine_mapping/simulation/susieinf/RunDirectory/R1_ref_n", n1, "_chr22_loci", i_loci, ".ld.gz")
  susieinf1 <- paste0(susieinf1, " --n ", n1)
  susieinf1 <- paste0(susieinf1, " --output-prefix /import/home/share/mingxuan/fine_mapping/simulation/susieinf/output/susieinf1_inflate_susieinf_loci_", i_loci, "_small_p", p, "_chr22_n", n1, "_corrected_loci", i_loci, "_", "_Ktrue", K_true, "_K",K,  "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                      "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i)
  system(susieinf1)


  cat(i, "-th rep finished.\n")
}





########################################################### pop1 ###########################################################
library(VCM)
library(XPASS)
library(mvtnorm)
library(susieR)
library(data.table)
library(RhpcBLASctl)
set.seed(1)
blas_set_num_threads(20)
source("/import/home/share/mingxuan/fine_mapping/utils.R")
source("/import/home/share/mingxuan/fine_mapping/XMAP.R")
source("/home/share/mingxuan/taam/code/sldxr.R")
bim <- fread("/home/share/mingxuan/wegene_qc2/height_ukb_qc2_chr22.bim")
bp_chr <- bim$V4

nrep <- 50

# n2_all <- c(10000, 15000, 20000)
n2 <- 20000
K_true <- 1#as.numeric(args[2]) #1
K <- 5

h1_omega <- 0.005
h2_omega <- 0.005
rho_omega <- 0.8

enrich <- 50
rho_sigma <- 0

pc_var1 <- 0.05
pc_var2 <- 0.2

loci <- fread(paste0("/home/share/mingxuan/fine_mapping/loci_chr22_1mb.txt"))
i_loci <- 29
idx_include <- (501:1000) + loci$idx_start[i_loci] - 1
p <- length(idx_include)

bim_loci <- bim[bim$V4 >= loci$left[i_loci] & bim$V4 <= loci$right[i_loci],]
bp_start <- bim_loci$V4[idx_include[1]]
bp_end <- bim_loci$V4[idx_include[p]]
snps <- bim_loci$V2[idx_include]



R2 <- data.matrix(fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_R2_ukb_in_sample_ref_n", n2, "_chr22_loci", i_loci, ".txt")))
R2 <- R2[idx_include, idx_include]
fwrite(R2, file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/susieinf/RunDirectory/R2_ref_n", n2, "_chr22_loci", i_loci, ".ld.gz"), col.names = F, quote = F, sep = "\t")
Beta2 <- data.matrix(fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_loci_", i_loci, "_small_p", p, "_Beta2_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                  "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".txt")))

for (i in 1:nrep) {

  z_scores2 <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_zs2_inflate_pc1_", pc_var2, "_corrected_loci_", i_loci, "_small_p", p, "_chr22_n", n2, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                            "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"))

  zs2 <- z_scores2$V2[loci$idx_left[i_loci]:loci$idx_right[i_loci]][idx_include]


  fwrite(data.frame(Z = zs2), file = paste0("/import/home/share/mingxuan/fine_mapping/simulation/susieinf/RunDirectory/simu_zs2_inflate_susieinf_loci_", i_loci, "_small_p", p, "_chr22_n", n2, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                            "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"))


  susieinf2 <- "python fine-mapping-inf/run_fine_mapping.py --z-col-name Z --save-tsv"
  susieinf2 <- paste0(susieinf2, " --sumstats /import/home/share/mingxuan/fine_mapping/simulation/susieinf/RunDirectory/simu_zs2_inflate_susieinf_loci_", i_loci, "_small_p", p, "_chr22_n", n2, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                      "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt")
  susieinf2 <- paste0(susieinf2, " --num-sparse-effects ", K)
  susieinf2 <- paste0(susieinf2, " --ld-file /import/home/share/mingxuan/fine_mapping/simulation/susieinf/RunDirectory/R2_ref_n", n2, "_chr22_loci", i_loci, ".ld.gz")
  susieinf2 <- paste0(susieinf2, " --n ", n2)
  susieinf2 <- paste0(susieinf2, " --output-prefix /import/home/share/mingxuan/fine_mapping/simulation/susieinf/output/susieinf2_inflate_susieinf_loci_", i_loci, "_small_p", p, "_chr22_n", n2, "_corrected_loci", i_loci, "_", "_Ktrue", K_true, "_K",K,  "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                      "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i)
  system(susieinf2)


  cat(i, "-th rep finished.\n")
}


