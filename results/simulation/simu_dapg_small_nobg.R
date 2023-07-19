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
n1 <- 20000
# n1 <- 10000
# n1 <- 15000
# n1 <- 20000
n2 <- 20000

K_all <- c(1:3)
# K_true <- 1
# K_true <- 2
# K_true <- 3


h1_omega <- 0.005
h2_omega <- 0.005
rho_omega <- 0.8

enrich <- 50
enrich2 <- enrich*10
rho_sigma <- 0

# per-snp variance for bg
omega1 <- h1_omega / p
omega2 <- h2_omega / p

# per-snp variance for causal
Sigma1 <- omega1 * enrich
Sigma2 <- omega2 * enrich

input_folder <- "/home/share/mingxuan/fine_mapping/simulation/dapg/RunDirectory/"

# LD matrix
R1 <- data.matrix(fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_R1_wg_in_sample_ref_n", n1, "_chr22_loci", i_loci, ".txt")))
R1_file <- paste0(
    input_folder,
    "dapg_input_nobg_loci_", i_loci, "_small_p", p, "_chr22_in_sample_n", n1, "_", n2,
    "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
    "_enrich", enrich2, "_rhoSigma", rho_sigma, "_nrep", nrep, ".LD1"
)
fwrite(R1[idx_include, idx_include],
    R1_file,
    sep = " ", col.names = F
)

# input file
z1_file <- paste0(
    input_folder,
    "dapg_input_nobg_loci_", i_loci, "_small_p", p, "_chr22_in_sample_n", n1, "_", n2,
    "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
    "_enrich", enrich2, "_rhoSigma", rho_sigma, "_nrep", nrep, ".z1"
)

for (K_true in K_all) {
    system(paste0("mkdir -p /home/share/mingxuan/fine_mapping/simulation/dapg/output/wg_n", n1, "_K", K_true, "_enrich",enrich2, "_nobg_default"))
    output_folder <- paste0("/home/share/mingxuan/fine_mapping/simulation/dapg/output/wg_n", n1, "_K", K_true, "_enrich",enrich2,"_nobg_default")

    for (i in 1:nrep) {
        z_scores1 <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_zs1_nobg_loci_",i_loci,"_small_p",p,"_chr22_n", n1, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                            "_enrich", enrich2, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"))
        z1 <- data.table(RSID = z_scores1[[1]], ZSCORE = z_scores1[[2]])
        write.table(z1, file = z1_file, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
        dapg <- paste0("./dap/dap_src/dap-g -d_z ", z1_file, " -d_ld ", R1_file, " -t 20 -o ", output_folder, "/dapg.", i)
        system(dapg)
        cat(i, "-th rep finished.\n")
    }
}














########################################################### pop2 ###########################################################
set.seed(1)
library(data.table)
library(RhpcBLASctl)
blas_set_num_threads(40)

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
# n1 <- 20000
n2 <- 20000

K_all <- c(1:3)
# K_true <- 1
# K_true <- 2
# K_true <- 3


h1_omega <- 0.005
h2_omega <- 0.005
rho_omega <- 0.8

enrich <- 50
enrich2 <- enrich*20
rho_sigma <- 0

# per-snp variance for bg
omega1 <- h1_omega / p
omega2 <- h2_omega / p

# per-snp variance for causal
Sigma1 <- omega1 * enrich
Sigma2 <- omega2 * enrich

input_folder <- "/home/share/mingxuan/fine_mapping/simulation/dapg/RunDirectory/"

# LD matrix
R2 <- data.matrix(fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_R2_ukb_in_sample_ref_n", n2, "_chr22_loci", i_loci, ".txt")))
R2_file <- paste0(
    input_folder,
    "dapg_input_nobg_loci_", i_loci, "_small_p", p, "_chr22_in_sample_n", n2,
    "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
    "_enrich2", enrich, "_rhoSigma", rho_sigma, "_nrep", nrep, ".LD2"
)
fwrite(R2[idx_include, idx_include],
    R2_file,
    sep = " ", col.names = F
)

# input file
z2_file <- paste0(
    input_folder,
    "dapg_input_nobg_loci_", i_loci, "_small_p", p, "_chr22_in_sample_n", n2,
    "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
    "_enrich", enrich2, "_rhoSigma", rho_sigma, "_nrep", nrep, ".z2"
)


for (K_true in K_all) {
    system(paste0("mkdir -p /home/share/mingxuan/fine_mapping/simulation/dapg/output/ukb_n", n2, "_K", K_true, "_enrich",enrich2, "_nobg_default"))
    output_folder <- paste0("/home/share/mingxuan/fine_mapping/simulation/dapg/output/ukb_n", n2, "_K", K_true, "_enrich",enrich2, "_nobg_default")

    for (i in 1:nrep) {
        z_scores2 <- fread(paste0("/import/home/share/mingxuan/fine_mapping/simulation/data/simu_zs2_nobg_loci_",i_loci,"_small_p",p,"_chr22_n", n2, "_loci", i_loci, "_", "_Ktrue", K_true, "_Omega", h1_omega, "_", h2_omega, "_", rho_omega,
                                "_enrich", enrich2, "_rhoSigma", rho_sigma, "_nrep", nrep, "_", i, ".txt"))

        z2 <- data.table(RSID = z_scores2[[1]], ZSCORE = z_scores2[[2]])
        write.table(z2, file = z2_file, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
        dapg <- paste0("./dap/dap_src/dap-g -d_z ", z2_file, " -d_ld ", R2_file, " -t 20 -o ", output_folder, "/dapg.", i)
        system(dapg)
        cat(i, "-th rep finished.\n")
    }
}



setwd("/home/share/mingxuan/fine_mapping/simulation/dapg/output/")
library(RhpcBLASctl)
blas_set_num_threads(40)

enrich <- 50
enrich2 <- enrich*10

n1_all <- c(5000, 10000, 15000, 20000)
n2 <- 20000
K_all <- c(1, 2, 3)

for (K in K_all) {
    for (i in c(1:50)) {
        # system(paste0('grep "((" ', "./ukb_n", n2, "_K", K, "_inflated_pc0.2/dapg.", i, " > ./ukb_n", n2, "_K", K, "_inflated_pc0.2/dapg_zs.", i))
        # system(paste0('grep "((" ', "./ukb_n", n2, "_K", K, "_msize_default/dapg.", i, " > ./ukb_n", n2, "_K", K, "_msize_default/dapg_zs.", i))
        system(paste0('grep "((" ', "./ukb_n", n2, "_K", K, "_enrich",enrich2, "_nobg_default/dapg.", i, " > ./ukb_n", n2, "_K", K, "_enrich",enrich2, "_nobg_default/dapg_zs.", i))
        # system(paste0('grep "((" ', "./ukb_n", n2, "_K", K, "/dapg.", i, " > ./ukb_n", n2, "_K", K, "/dapg_zs.", i))
        for (n1 in n1_all) {
            # system(paste0('grep "((" ', "./wg_n", n1, "_K", K, "_inflated_pc0.05/dapg.", i, " > ./wg_n", n1, "_K", K, "_inflated_pc0.05/dapg_zs.", i))
            # system(paste0('grep "((" ', "./wg_n", n1, "_K", K, "_msize_default/dapg.", i, " > ./wg_n", n1, "_K", K, "_msize_default/dapg_zs.", i))
            system(paste0('grep "((" ', "./wg_n", n1, "_K", K, "_enrich",enrich2, "_nobg_default/dapg.", i, " > ./wg_n", n1, "_K", K, "_enrich",enrich2, "_nobg_default/dapg_zs.", i))
            # system(paste0('grep "((" ', "./wg_n", n1, "_K", K, "/dapg.", i, " > ./wg_n", n1, "_K", K, "/dapg_zs.", i))
        }
    }
}