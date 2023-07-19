#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
########################################################### pop1 ###########################################################
library(mvtnorm)
library(data.table)
set.seed(1)
source("/import/home/share/mingxuan/fine_mapping/utils.R")
source("/import/home/share/mingxuan/fine_mapping/XMAP.R")
source("/home/share/mingxuan/taam/code/sldxr.R")
# system(paste0("export MKL_NUM_THREADS=",nthreads))
# system(paste0("export NUMEXPR_NUM_THREADS=",nthreads))
# system(paste0("export OMP_NUM_THREADS=",nthreads))
# export MKL_NUM_THREADS=10
# export NUMEXPR_NUM_THREADS=10
# export OMP_NUM_THREADS=10

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

pc_var1 <- as.numeric(args[4])#0.2
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

LD1_file <- paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/paintor/working_input/paintor_mcmc_input_inflate_pc1_", pc_var1, "_", pc_var2, "_loci_", i_loci, "_small_p", p, "_chr22_paintor_in_sample_n", n1, "_", n2,
                   "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                   "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".LD1")
LD2_file <- paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/paintor/working_input/paintor_mcmc_input_inflate_pc1_", pc_var1, "_", pc_var2, "_loci_", i_loci, "_small_p", p, "_chr22_paintor_in_sample_n", n1, "_", n2,
                   "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                   "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".LD2")
LD_input_file <- paste0("/home/share/mingxuan/fine_mapping/simulation/data/mscaviar/input/mscaviar_LD_input_inflate_pc1_", pc_var1, "_", pc_var2, "_loci_", i_loci, "_small_p", p, "_chr22_paintor_in_sample_n", n1, "_", n2,
                        "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                        "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep)
write.table(rbind(LD1_file, LD2_file),
            file = LD_input_file,
            col.names = F, row.names = F, sep = " ", quote = F)

ofile <- paste0("/home/share/mingxuan/fine_mapping/simulation/data/mscaviar/output/mscaviar_output_inflate_pc1_", pc_var1,"_",pc_var2, "_loci_", i_loci, "_small_p", p, "_chr22_paintor_in_sample_n", n1, "_", n2,
                "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", 1:nrep)
z1_file <- paste0("/home/share/mingxuan/fine_mapping/simulation/data/mscaviar/input/mscaviar_Zs1_inflate_pc1_", pc_var1, "_loci_", i_loci, "_small_p", p, "_chr22_paintor_in_sample_n", n1, "_", n2,
                  "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                  "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", 1:nrep)
z2_file <- paste0("/home/share/mingxuan/fine_mapping/simulation/data/mscaviar/input/mscaviar_Zs2_inflate_pc1_", pc_var2, "_loci_", i_loci, "_small_p", p, "_chr22_paintor_in_sample_n", n1, "_", n2,
                  "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                  "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", 1:nrep)
z_input_file <- paste0("/home/share/mingxuan/fine_mapping/simulation/data/mscaviar/input/mscaviar_Z_input_inflate_pc1_", pc_var1,"_",pc_var2, "_loci_", i_loci, "_small_p", p, "_chr22_paintor_in_sample_n", n1, "_", n2,
                       "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                       "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", 1:nrep)


for (i in 1:nrep) {
  z_scores1 <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_zs1_inflate_pc1_", pc_var1, "_loci_", i_loci, "_small_p", p, "_chr22_n", n1, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                            "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"))
  z_scores2 <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_zs2_inflate_pc1_", pc_var2, "_loci_", i_loci, "_small_p", p, "_chr22_n", n2, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                            "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"))

  zs1 <- z_scores1[loci$idx_left[i_loci]:loci$idx_right[i_loci],][idx_include,]
  zs2 <- z_scores2[loci$idx_left[i_loci]:loci$idx_right[i_loci],][idx_include,]

  write.table(zs1,
              file = z1_file[i],
              col.names = F, row.names = F, sep = " ", quote = F)
  write.table(zs2,
              file = z2_file[i],
              col.names = F, row.names = F, sep = " ", quote = F)
  write.table(rbind(z1_file[i], z2_file[i]),
              file = z_input_file[i],
              col.names = F, row.names = F, sep = " ", quote = F)

  mscaviar <- paste0("/home/mcaiad/MsCAVIAR/MsCAVIAR -l ", LD_input_file, " -z ", z_input_file[i], " -n ", n1, ",", n2, " -c ",K_true," -o ", ofile[i])
  system(mscaviar)
  cat(i, "-th rep finished.\n")
}

