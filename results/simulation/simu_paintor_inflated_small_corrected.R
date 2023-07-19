#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
########################################################### pop1 ###########################################################
library(XPASS)
library(mvtnorm)
library(susieR)
library(data.table)
library(RhpcBLASctl)
set.seed(1)
blas_set_num_threads(10)
source("/import/home/share/mingxuan/fine_mapping/utils.R")
source("/import/home/share/mingxuan/fine_mapping/XMAP.R")
source("/home/share/mingxuan/taam/code/sldxr.R")
bim <- fread("/home/share/mingxuan/wegene_qc2/height_ukb_qc2_chr22.bim")
bp_chr <- bim$V4

nrep <- 50

# n1_all <- c(5000, 10000, 15000, 20000)
# n2_all <- c(20000)
n1 <- as.numeric(args[1]) #5000
n2 <- as.numeric(args[2])#20000
K_true <- as.numeric(args[3])
K <- K_true

h1_omega <- 0.005
h2_omega <- 0.005
rho_omega <- 0.8

enrich <- 50
rho_sigma <- 0

pc_var1 <- as.numeric(args[4]) #0.2
pc_var2 <- as.numeric(args[5])#0

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
fwrite(R1[idx_include, idx_include], paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/paintor/working_input/paintor_mcmc_input_inflate_pc1_", pc_var1, "_", pc_var2, "_corrected_loci_", i_loci, "_small_p", p, "_chr22_paintor_in_sample_n", n1, "_", n2,
                                            "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                            "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".LD1"),
       sep = " ", col.names = F)
R2 <- data.matrix(fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_R2_ukb_in_sample_ref_n", n2, "_chr22_loci", i_loci, ".txt")))
fwrite(R2[idx_include, idx_include], paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/paintor/working_input/paintor_mcmc_input_inflate_pc1_", pc_var1, "_", pc_var2, "_corrected_loci_", i_loci, "_small_p", p, "_chr22_paintor_in_sample_n", n1, "_", n2,
                                            "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                            "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".LD2"),
       sep = " ", col.names = F)
# system(paste0("cp /import/home/share/mingxuan/fine_mapping/simulation/data/simu_R1_wg_in_sample_ref_n", n1, "_chr22_loci", i_loci, ".txt /import/home/share/mingxuan/fine_mapping/simulation/data/paintor/working_input/paintor_mcmc_input_loci_34_chr22_paintor_in_sample_n", n1, "_", n2,
#               "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
#               "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".LD1"))
# system(paste0("cp /import/home/share/mingxuan/fine_mapping/simulation/data/simu_R2_ukb_in_sample_ref_n", n2, "_chr22_loci", i_loci, ".txt /import/home/share/mingxuan/fine_mapping/simulation/data/paintor/working_input/paintor_mcmc_input_loci_34_chr22_paintor_in_sample_n", n1, "_", n2,
#               "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
#               "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".LD2"))

annot <- paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/paintor/working_input/paintor_mcmc_input_inflate_pc1_", pc_var1, "_", pc_var2, "_corrected_loci_", i_loci, "_small_p", p, "_chr22_paintor_in_sample_n", n1, "_", n2,
                "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".annotations")
dummy <- data.table(dummy = rep(1, p))
fwrite(dummy, file = annot, sep = " ", col.names = TRUE)

input_folder <- "/import/home/share/mingxuan/fine_mapping/simulation/data/paintor/working_input/"
zfile <- paste0("paintor_mcmc_input_inflate_pc1_", pc_var1, "_", pc_var2, "_corrected_loci_", i_loci, "_small_p", p, "_chr22_paintor_in_sample_n", n1, "_", n2,
                "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep)
input_file <- paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/paintor/working_input/files_apaintor_mcmc_input_inflate_pc1_", pc_var1, "_", pc_var2, "_corrected_loci_", i_loci, "_small_p", p, "_chr22_paintor_in_sample_n", n1, "_", n2,
                     "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                     "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep)
for (i in 1:nrep) {
  z_scores1 <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_zs1_inflate_pc1_", pc_var1, "_corrected_loci_", i_loci, "_small_p", p, "_chr22_n", n1, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                            "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"))
  z_scores2 <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_zs2_inflate_pc1_", pc_var2, "_corrected_loci_", i_loci, "_small_p", p, "_chr22_n", n2, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                            "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"))

  zs1 <- z_scores1[loci$idx_left[i_loci]:loci$idx_right[i_loci],][idx_include,]
  zs2 <- z_scores2[loci$idx_left[i_loci]:loci$idx_right[i_loci],][idx_include,]

  z_paintor <- data.table(RSID = zs1[[1]], ZSCORE.P1 = zs1[[2]], ZSCORE.P2 = zs2[[2]])
  fwrite(z_paintor, paste0(input_folder, zfile), sep = " ", col.names = TRUE)


  write.table(zfile, file = input_file, sep = " ", col.names = FALSE, row.names = F, quote = F)

  paintor <- paste0("/home/mcaiad/PAINTOR_V3.0/PAINTOR -input ", input_file, " -in ", input_folder, " -out /import/home/share/mingxuan/fine_mapping/simulation/data/paintor/output/ -Zhead ZSCORE.P1,ZSCORE.P2 -LDname LD1,LD2 -enumerate ", K, " -annotations dummy -RESname ", i, "_exact_results")
  system(paintor)
  cat(i, "-th rep finished.\n")
}

