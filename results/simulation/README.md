# Simulation in XMAP paper
Here, we provide the codes for generating the simulation results in XMAP paper.
The file names are specified with their functionalities + specific simulation scenario in XMAP paper.

## Simulation scenario
- Polygenic + causal effects from normal distributions: *_small_chr22.R
- Polygenic + causal effects from scaled t distributions: *_small_t_chr22.R
- Causal effects + counfounding (introduced by PC): *_inflated_small.R
- Causal effects only: *_small_nobg.R

## Codes for generating simulation data
- Codes for generating effect sizes: Start with 'construct_effect_size'
- Codes for generating GWAS data: Start with 'construct_gwas'

## Codes for running fine-mapping methods
- XMAP: Start with "simu_XMAP"
- SuSiEx: Start with "simu_susie"
- MsCAVIAR: Start with "simu_mscaviar"
- PAINTOR: Start with "simu_paintor"
- SuSiE-inf: Start with "simu_susieinf"
- SuSiE: Start with "simu_susie"
- FINEMAP: Start with "simu_fm"
- DAP-G: Start with "simu_dapg"

