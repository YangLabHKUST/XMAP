# XMAP
[![DOI](https://zenodo.org/badge/598050665.svg)](https://zenodo.org/badge/latestdoi/598050665)

XMAP is a computationally efficient and statistically accurate method for fine-mapping causal variants using GWAS summary statistics. 

Check out our manuscript in Nature Communications:

- [Nature Communications](https://doi.org/10.1038/s41467-023-42614-7)
- [Full text](https://urldefense.com/v3/__https://rdcu.be/dpGfM__;!!KjDnqvtInNPT!gO0tXuDsQM0696a74hDXc9wMM_F1ke7-AswOo6Wj32PUCbvTRjP5nd_lTRKx-GDr66NJcCinlJEa6cn1BzyE9Cz2tasbsgWzSCaT8g$)
- [bioRxiv preprint](https://doi.org/10.1101/2023.03.30.534832)

Briefly, it can leverage different LD structures of genetically diverged populations to better distinguish causal variants from a set of associated variants. By jointly modeling SNPs with putative causal effects and polygenic effects, XMAP allows a linear-time computational cost to identify multiple causal variants, even in the presence of an over-specified number of causal variants. It further corrects confounding bias hidden in the GWAS summary statistics to reduce false positive findings and improve replication rates.

The fine-mapping results given by XMAP can be further used for downstream analysis to illuminate the causal mechanisms at different cascades of biological processes, including tissues, cell populations, and individual cells. In particular, XMAP results can be effectively integrated with single-cell datasets to identify disease/trait-relevant cells.
![XMAP_overview](https://github.com/YangLabHKUST/XMAP/blob/main/results/flowchart.png)

# Installation

* Prerequisites: XMAP is developed under R (version >= 3.6.1).

* Latest version: The latest developmental version of XMAP can be downloaded from GitHub and installed from source by 
```
# Install devtools, if necessary
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
# Install XMAP
devtools::install_github("YangLabHKUST/XMAP")
# load XMAP
library(XMAP)
```

# Manual
Please see https://mxcai.github.io/XMAP-tutorial/index.html for details. In the R terminal, please use the command `?XMAP` to access the help documents.

# Tutorial
We provide a [tutorial website](https://mxcai.github.io/XMAP-tutorial/index.html) of XMAP for analyzing cross-population GWAS data. Please see the vignettes in the tutorial for details of installation, usage and visualization. The datasets involved in the tutorial can be downloaded [here](https://hkustconnect-my.sharepoint.com/:f:/g/personal/mcaiad_connect_ust_hk/EhJHXBkK_DNBjLFbIPjMeaoBFlmFwlz0F_uXXU0kvIrVGg?e=sTEh8O).

# Reproducibility
We provide codes and fine-mapping results presented in the XMAP manuscript [here](https://github.com/YangLabHKUST/XMAP/tree/main/results)

# Operating systems tested on:
macOS Ventura 13.0 


Windows 10 Enterprise Version


Ubuntu 18.04.5 LTS (Bionic Beaver)

# License
XMAP is licensed under the GNU General Public License v3.0.

# Reference
Mingxuan Cai, Zhiwei Wang, Jiashun Xiao, Xianghong Hu, Gang Chen, Can Yang (2023). XMAP: Cross-population fine-mapping by leveraging genetic diversity and accounting for confounding bias. Nature Communications, 14(1),6870. [Link](https://doi.org/10.1038/s41467-023-42614-7)

# Contact
Improvements and new features of XMAP will be updated on a regular basis. Please feel free to contact Mingxuan Cai (mx.cai@cityu.edu.hk) or Prof. Can Yang (macyang@ust.hk) if any questions. 
