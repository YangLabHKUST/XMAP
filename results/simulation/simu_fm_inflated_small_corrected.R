########################################################### pop1 ###########################################################
set.seed(1)
library(data.table)
library(RhpcBLASctl)
blas_set_num_threads(10)

nrep <- 50

bim <- fread("/home/share/mingxuan/wegene_qc2/height_ukb_qc2_chr22.bim")
bp_chr <- bim$V4

loci <- fread(paste0("/home/share/mingxuan/fine_mapping/loci_chr22_1mb.txt"))
i_loci <- 29
idx_include <- (501:1000) + loci$idx_start[i_loci] - 1
p <- length(idx_include)

# n1_all <- c(5000, 10000, 15000, 20000)
# n2_all <- c(20000)
# n1 <- 5000
# n1 <- 10000
# n1 <- 15000
n1 <- 20000

K_all <- c(1:3)
# K_true <- 2
# K_true <- 3


h1_omega <- 0.005
h2_omega <- 0.005
rho_omega <- 0.8

enrich <- 50
rho_sigma <- 0
pc_var1 <- 0.05

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

input_folder <- "/home/share/mingxuan/fine_mapping/simulation/finemap/RunDirectory/"

# LD matrix
R1 <- data.matrix(fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_R1_wg_in_sample_ref_n", n1, "_chr22_loci", i_loci, ".txt")))
R1_file <- paste0(
  input_folder,
  "finemap_input_inflate_pc1_", pc_var1, "_corrected_loci_", i_loci, "_small_p", p, "_chr22_in_sample_n", n1,
  "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
  "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".ld"
)
fwrite(R1[idx_include, idx_include],
       R1_file,
       sep = " ", col.names = F
)

# input file
z1_file <- paste0(
  input_folder,
  "finemap_input_inflate_pc1_", pc_var1, "_corrected_loci_", i_loci, "_small_p", p, "_chr22_in_sample_n", n1,
  "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
  "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".z"
)


for (K_true in K_all) {
  master_file <- paste0(
    input_folder,
    "wg_n",
    n1, "_K", K_true,
    ".master"
  )

  system(paste0("mkdir -p /home/share/mingxuan/fine_mapping/simulation/finemap/output/wg_n", n1, "_K", K_true, "_enrich", enrich, "_inflated_pc", pc_var1, "_corrected_default"))
  output_folder <- paste0("/home/share/mingxuan/fine_mapping/simulation/finemap/output/wg_n", n1, "_K", K_true, "_enrich", enrich, "_inflated_pc", pc_var1, "_corrected_default/")

  for (i in 1:nrep) {
    z_scores1 <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_zs1_inflate_pc1_", pc_var1, "_corrected_loci_", i_loci, "_small_p", p, "_chr22_n", n1, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                              "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"))
    se <- 1 / sqrt(n1)
    z1 <- data.table(rsid = z_scores1[[1]][loci$idx_left[i_loci]:loci$idx_right[i_loci]][idx_include], chromosome = rep(22, p), position = c(1:p), allele1 = rep("A", p), allele2 = rep("T", p), maf = rep(0.05, p), beta = z_scores1[[2]][loci$idx_left[i_loci]:loci$idx_right[i_loci]][idx_include] * se, se = rep(se, p))
    write.table(z1, file = z1_file, sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)

    ms <- data.frame(z = z1_file, ld = R1_file, snp = paste0(output_folder, i, ".snp"), config = paste0(output_folder, i, ".config"), cred = paste0(output_folder, i, ".cred"), log = paste0(output_folder, i, ".log"), n_samples = n1)
    write.table(ms, file = master_file, sep = ";", col.names = TRUE, row.names = FALSE, quote = FALSE)

    fm <- paste0("./finemap_v1.4_x86_64/finemap_v1.4_x86_64 --sss --in-files ", master_file, " --n-threads 20")
    system(fm)
    cat(i, "-th rep finished.\n")
  }
}


########################################################### pop2 ###########################################################
set.seed(1)
library(data.table)
library(RhpcBLASctl)
blas_set_num_threads(20)

nrep <- 50

bim <- fread("/home/share/mingxuan/wegene_qc2/height_ukb_qc2_chr22.bim")
bp_chr <- bim$V4

loci <- fread(paste0("/home/share/mingxuan/fine_mapping/loci_chr22_1mb.txt"))
i_loci <- 29
idx_include <- (501:1000) + loci$idx_start[i_loci] - 1
p <- length(idx_include)

n2 <- 20000

K_all <- c(1:3)
# K_true <- 2
# K_true <- 3


h1_omega <- 0.005
h2_omega <- 0.005
rho_omega <- 0.8
pc_var2 <- 0.2

loci <- fread(paste0("/home/share/mingxuan/fine_mapping/loci_chr22_1mb.txt"))
i_loci <- 29
idx_include <- (501:1000) + loci$idx_start[i_loci] - 1
p <- length(idx_include)

enrich <- 50
rho_sigma <- 0

# per-snp variance for bg
omega1 <- h1_omega / p
omega2 <- h2_omega / p

# per-snp variance for causal
Sigma1 <- omega1 * enrich
Sigma2 <- omega2 * enrich

input_folder <- "/home/share/mingxuan/fine_mapping/simulation/finemap/RunDirectory/"


# LD matrix
R2 <- data.matrix(fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_R2_ukb_in_sample_ref_n", n2, "_chr22_loci", i_loci, ".txt")))
R2_file <- paste0(
  input_folder,
  "finemap_input_inflate_pc1_", pc_var2, "_corrected_loci_", i_loci, "_small_p", p, "_chr22_in_sample_n", n2,
  "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
  "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".ld"
)


fwrite(R2[idx_include, idx_include],
       R2_file,
       sep = " ", col.names = F
)


# input file
z2_file <- paste0(
  input_folder,
  "finemap_input_inflate_pc1_", pc_var2, "_corrected_loci_", i_loci, "_small_p", p, "_chr22_in_sample_n", n2,
  "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
  "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".z"
)


for (K_true in K_all) {
  master_file <- paste0(
    input_folder,
    "ukb_n",
    n2, "_K", K_true,
    ".master"
  )

  system(paste0("mkdir -p /home/share/mingxuan/fine_mapping/simulation/finemap/output/ukb_n", n2, "_K", K_true, "_enrich", enrich, "_inflated_pc", pc_var2, "_corrected_default"))
  output_folder <- paste0("/home/share/mingxuan/fine_mapping/simulation/finemap/output/ukb_n", n2, "_K", K_true, "_enrich", enrich, "_inflated_pc", pc_var2, "_corrected_default/")

  for (i in 1:nrep) {

    z_scores2 <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_zs2_inflate_pc1_", pc_var2, "_corrected_loci_", i_loci, "_small_p", p, "_chr22_n", n2, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                              "_enrich", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"))
    se <- 1 / sqrt(n2)
    z2 <- data.table(rsid = z_scores2[[1]][loci$idx_left[i_loci]:loci$idx_right[i_loci]][idx_include], chromosome = rep(29, p), position = c(1:p), allele1 = rep("A", p), allele2 = rep("T", p), maf = rep(0.05, p), beta = z_scores2[[2]][loci$idx_left[i_loci]:loci$idx_right[i_loci]][idx_include] * se, se = rep(se, p))
    write.table(z2, file = z2_file, sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)

    ms <- data.frame(z = z2_file, ld = R2_file, snp = paste0(output_folder, i, ".snp"), config = paste0(output_folder, i, ".config"), cred = paste0(output_folder, i, ".cred"), log = paste0(output_folder, i, ".log"), n_samples = n2)
    write.table(ms, file = master_file, sep = ";", col.names = TRUE, row.names = FALSE, quote = FALSE)

    fm <- paste0("./finemap_v1.4_x86_64/finemap_v1.4_x86_64 --sss --in-files ", master_file, " --n-threads 20")
    system(fm)
    cat(i, "-th rep finished.\n")
  }
}