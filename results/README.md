# Results in XMAP paper
Here, we provide the analysis results and codes for generating them in XMAP paper.

## Codes
### Figure generating codes in html:
Figure 4: [html](https://github.com/YangLabHKUST/XMAP/blob/main/results/Fig4_LDL_ncausal.html); [notebook](https://github.com/YangLabHKUST/XMAP/blob/main/results/Fig4_LDL_ncausal.ipynb)

Figure 5: [html](https://github.com/YangLabHKUST/XMAP/blob/main/results/Fig5_replicate.html); [notebook](https://github.com/YangLabHKUST/XMAP/blob/main/results/Fig5_replicate.ipynb)

Figure 6: [html](https://github.com/YangLabHKUST/XMAP/blob/main/results/Fig6_height_multiple_signal.html); [notebook](https://github.com/YangLabHKUST/XMAP/blob/main/results/Fig6_height_multiple_signal.ipynb)

Figure 7: [html](https://github.com/YangLabHKUST/XMAP/blob/main/results/Fig7_SCAVENGE.ipynb); [notebook](https://github.com/YangLabHKUST/XMAP/blob/main/results/Fig7_SCAVENGE.ipynb)

### LDL
Fine-mapping for all loci on the genome: LDL_all_3pop.R

Extract PIP from the fine-mapping output: pip_LDL.R

### height
Fine-mapping for all loci on the genome: height_all_data.R

Extract PIP from the fine-mapping output: pip_height_all_data.R

### 12 blood traits
Fine-mapping for all loci on the genome: All_blood_traits_BBJ_UKB.R

Extract PIP from the fine-mapping output: pip_blood.R

SCAVENGE analysis for identifying trait-relevant cells: scavenge_analysis.R

Generate Figure 7 in the main text: plot_scavenge_blood.R

## Fine-mapping Results
### LDL
Link: https://hkustconnect-my.sharepoint.com/:f:/g/personal/mcaiad_connect_ust_hk/ElvWk3EXEEVCpej0aJM2dx0B0mECe8b_L18KPdinu8UXEA?e=09G3ZE

XMAP (EUR+AFR+EAS): LDL_3pop_PIP.txt

XMAP (EUR+AFR): LDL_EUR_AFR_PIP.txt

XMAP (EUR+EAS): LDL_EUR_EAS_PIP.txt

XMAP (AFR+EAS): LDL_AFR_EAS_PIP.txt

SuSiE (EUR): LDL_EUR_PIP.txt

SuSiE (AFR): LDL_AFR_PIP.txt

SuSiE (EAS): LDL_EAS_PIP.txt

### height
Link: https://hkustconnect-my.sharepoint.com/:f:/g/personal/mcaiad_connect_ust_hk/Eptgw2McNLRKjsypetdRwzABaA_JM20jEeGgSwwc_GfK1g?e=XNloZD
#### Discovery cohorts:
XMAP (UKBB+Chinese): height_UKB_WG_K10_PIP.txt  height_UKB_WG_K15_PIP.txt

XMAP with C=I (UKBB+Chinese): height_UKB_WG_CI_K10_PIP.txt  height_UKB_WG_CI_K15_PIP.txt

SuSiE (UKBB): height_UKB_K10_PIP.txt  height_UKB_K15_PIP.txt

SuSiE (Chinese): height_WG_K10_PIP.txt  height_WG_K15_PIP.txt

#### Replication cohorts:
XMAP (Sibship+Chinese): height_Sibship_WG_PIP.txt

SuSiE (Sibship): height_Sibship_PIP.txt

SuSiE (BBJ): height_BBJ_PIP.txt

### 12 blood traits
Link: https://hkustconnect-my.sharepoint.com/:f:/g/personal/mcaiad_connect_ust_hk/EiexnAhCaEZJj-1A3EhOpwQBAJZPF6eMYuiMQf0QH8flzg?e=NErpj7

XMAP (UKBB+BBJ): *_ukb_bbj_PIP.txt

SuSiE (UKBB): *_ukb_PIP.txt

SuSiE (BBJ): *_bbj_PIP.txt

SCAVENGE analysis of blood traits: https://hkustconnect-my.sharepoint.com/:f:/g/personal/mcaiad_connect_ust_hk/EkcuIxsmluFDrFcWj7C3Ly8BfrWwJ-LeuQZxC-MODjbXUw?e=Rh7ZPv

where * = WBC, Lym, Neutro, Mono, Eosino, Baso, Plt, RBC, MCH, MCHC, MCV, Hb
